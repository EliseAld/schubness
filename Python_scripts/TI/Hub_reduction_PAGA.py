from skhubness.neighbors import kneighbors_graph
import scanpy as sc
import anndata
import os

# load data
def load_preprocess(path, n_dim):
    adata = anndata.read_h5ad(path)
    # reduce the dimension
    sc.pp.recipe_zheng17(adata)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_dim)
    print("We have "+str(adata.obsm['X_pca'].shape[0])+" cells over "+str(adata.obsm['X_pca'].shape[1])+" PCs")
    return(adata)

def prune_dict(dictionary, mask):
    result = {}
    for key in dictionary:
        if isinstance(dictionary[key], dict):
            newmask = [maskname[mask.find(".") + 1:] for maskname in mask if maskname.startswith(key + ".")]
            result[k] = filterstage1(dictionary[key], newmask)
        elif key in mask:
            result[key] = dictionary[key]
    return result

# Apply hub reduction and export the graph
def hub_paga(adata, methods_set,path_out):
    hubness_choice = {'nothing':(None,None),
                   'mp_normal':('mp',{'method': 'normal'}),
                   'ls':('ls',None),
                   'ls_nicdm':('ls',{'method': 'nicdm'}),
                   'dsl':('dsl',None)
                   }
    hubness_methods = prune_dict(hubness_choice,methods_set)
    all_adata = dict()
    # get the different graphs
    for method_name, (hubness, hubness_params) in hubness_methods.items():
        all_adata[method_name] = adata.copy()
        all_adata[method_name].obsm[method_name] = kneighbors_graph(all_adata.get(method_name).obsm['X_pca'],
                                                                n_neighbors=10,
                                                                hubness=hubness,
                                                                hubness_params=hubness_params)
    # Apply PAGA on the different graphs
    for method_name in all_adata.keys():
        sc.pp.neighbors(all_adata[method_name], n_neighbors=10, use_rep=method_name)
        sc.tl.leiden(all_adata[method_name])
        sc.tl.paga(all_adata[method_name], groups="leiden")
        #sc.pl.paga(all_adata.get(method_name), color='leiden')
        #print(method_name+' completed')
    # remove graphs
    for method_name in all_adata.keys():
        del(all_adata[method_name].obsm[method_name])
    # Export PAGA results
    for method_name in all_adata.keys():
        all_adata[method_name].write_h5ad(filename=path_out+method_name+".h5ad")
    print("export completed")

if __name__ == '__main__':

    if input('Do you have the h5ad already? y/n \n').lower() == 'n':
        main_path = "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/"
        print('Available data sets:')
        print([f.name for f in os.scandir(main_path) if f.is_dir()])
        dataset_folder = input("Which dataset do you want? \n")
        print('Available expression matrices:')
        print([f for f in os.listdir(main_path+dataset_folder) if f.endswith(".h5ad")])
        path = main_path+dataset_folder+'/data.h5ad'
        n_dim = int(input("How many PCs do you mant? \n"))
        adata = load_preprocess(path, n_dim)
        methods_choice = ['nothing','mp_normal','ls','ls_nicdm','dsl']
        print("Hub reduction methods are "+str(methods_choice))
        choice1 = input("Do you mant to test all methods? y/n\n")
        if choice1=='y':
            methods_set = methods_choice
        else:
            methods_set = input("Which methods? \n")
        path_out = main_path+dataset_folder+'/'+str(n_dim)+'_dims_'
        hub_paga(adata,methods_set,path_out)
    else:
        print("You can go to R directly")
