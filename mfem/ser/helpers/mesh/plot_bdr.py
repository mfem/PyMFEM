import numpy as np
import ifigure.utils.geom
from ifigure.interactive import figure, solid
import matplotlib.cm as cm
from .mesh_viewer import open_meshviewer

def plot_bdr_mesh(meshname = '', mesh = None, idx = 'all',
                  viewer = None, cmap='Dark2', name = '', return_data = False):
    if mesh is None:  mesh = globals()[meshname]
    if viewer is None:  viewer = open_meshviewer()
    
    nbdr = mesh.bdr_attributes.Size()
    kbdr = mesh.bdr_attributes.ToList()
    nbe = mesh.GetNBE()

    
    idxarr = kbdr if str(idx) == 'all' else idx    
    attr = np.array([mesh.GetBdrElement(i).GetAttribute() for i in range(nbe)])
    ivert = np.vstack([mesh.GetBdrElement(i).GetVerticesArray() for i in range(nbe)])
#    print [x for x in ivert[attr == 1, :]]

    cmap = cm.get_cmap(cmap)
    for ii in idxarr:
        data = np.dstack([np.vstack([mesh.GetVertexArray(int(k)) for k in x]) 
                      for x in ivert[attr == ii, :]])
        data = np.rollaxis(data, 2, 0)  ### [i_triangle, i_vert, xyz]    
        #data[:,:,2] =         data[:,:,2]*100
        obj= viewer.solid(data, facecolor = cmap(float(ii)/nbdr),
                          edgecolor=(0,0,0,1))
        if obj is not None: obj.rename(name + 'bdry_'+str(ii))
    if return_data: return data

def plot_bdr(meshname = '', mesh = None, idx = 'all', viewer = None,
             name = '', return_data = False,  **kwargs):
    if mesh is None:  mesh = globals()[meshname]
    if viewer is None:  viewer = open_meshviewer()

    nbdr = mesh.bdr_attributes.Size()
    kbdr = mesh.bdr_attributes.ToList()
    nbe = mesh.GetNBE()

    idxarr = kbdr if str(idx) == 'all' else idx
    attr = np.array([mesh.GetBdrElement(i).GetAttribute() for i in range(nbe)])
#    ivert = np.vstack([mesh.GetBdrElement(i).GetVerticesArray() for i in range(nbe)])
#    print [x for x in ivert[attr == 1, :]]

    for ii in idxarr:
        edges = np.array([mesh.GetBdrElementEdges(i)[0] for i in range(nbe)
                          if mesh.GetBdrElement(i).GetAttribute()==ii]).flatten()
        d = {}
        for x in edges:d[x] = d.has_key(x)
        edges = [x for x in d.keys() if not d[x]]
        ivert = [mesh.GetEdgeVertices(x) for x in edges]
        ivert = ifigure.utils.geom.connect_pairs(ivert)
        vv = np.vstack([mesh.GetVertexArray(i) for i in ivert])
        obj = viewer.plot(vv[:,0], vv[:,1], vv[:,2], **kwargs)
        if obj is not None: obj.rename(name + 'bdry_'+str(ii))        

    if return_data: return vv

def plot_domain(meshname = '', mesh = None, idx = 'all',
                viewer = None, cmap='Dark2', return_data = False):
    if mesh is None:  mesh = globals()[meshname]
    if viewer is None:  viewer = open_meshviewer()
    cmap = cm.get_cmap(cmap)    

    nbdr = mesh.bdr_attributes.Size()
    ndom = mesh.attributes.Size()    
    kdom = mesh.attributes.ToList()
    nbe = mesh.GetNBE()
    nel = mesh.GetNE()    

    idxarr = kdom if str(idx) == 'all' else idx    
    dom_attr = np.array([mesh.GetElement(i).GetAttribute() for i in range(nel)])
    bdr_attr = np.array([mesh.GetBdrElement(i).GetAttribute() for i in range(nbe)])
    dfaces = np.array([mesh.GetElementFaces(i)[0] for i in range(nel)])
    bfaces = [mesh.GetBdrElementFace(i)[0] for i in range(nbe)]
    bfaces0 = np.unique(bfaces)
    viewer.update(False)
    for ii in idxarr:
        domain_faces = np.unique(np.array([dfaces[dom_attr==ii]]).flatten())
        faces = np.intersect1d(domain_faces, bfaces0)
        bindex = [bfaces.index(f) for f in faces]
        bindex = np.unique(bdr_attr[bindex])
        
        data = np.dstack([np.vstack([mesh.GetVertexArray(x)
                                     for x in mesh.GetFaceVertices(kk)]) for kk in faces])
        data = np.rollaxis(data, 2, 0)  ### [i_triangle, i_vert, xyz]    
        #data[:,:,2] =         data[:,:,2]*100
        col =cmap(float(ii)/ndom)
        obj= viewer.solid(data, facecolor = col, 
                          linewidth = 0.0)
        if obj is not None: obj.rename('domain_'+str(ii))                
        plot_bdr(mesh = mesh, idx = bindex, viewer = viewer,  color = [0,0,0,1],
                 name = 'domain_'+str(ii)+'_')
    viewer.update(True)        
    if return_data: return data

    
