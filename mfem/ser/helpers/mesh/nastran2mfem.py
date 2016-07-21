'''
    read .nas file and make MFEM mesh file

    note: it reads only minimum set of grid data

    supported card
        GRID
        CTETRA
        CHEXA
        CTRIA6
        CQUAD8
        PSOLID
        PSHELL

    index in nas starts from 1. 
    index in mfem starts from 0.
'''
import numpy as np

class NASReader(object):
    def __init__(self, filename):
        self.filename = filename
        self.dataset = None
    def load(self):
        fid = open(self.filename, 'r')

        cl = ''

        while True:
            l = self._read_line(fid)
            if l.startswith('$ Grid data section'): break

        grids = []  ### reading grid      
        while True:
            l = self._read_line(fid)
            if l.startswith('$ Element data section'):break
            grids.append(self.parse_grid_fixed(l))

        grids = [np.array((float(g[3]), float(g[4]), float(g[5]))) for g in grids]
        grids = np.vstack(grids)

        elems = {'TRIA6':[],
                'TETRA':[],
                'HEXA':[],
                'QUAD8':[]} ### reading elements
        while True:            
            l = self._read_line(fid)
            if l.startswith('$ Property data section'):break
            if l.startswith('CTRIA6'):
                elems['TRIA6'].append(self.parse_tria6_fixed(l))
            elif l.startswith('CTETRA'):
                elems['TETRA'].append(self.parse_tetra_fixed(l))
            elif l.startswith('CHEXA'):
                elems['HEXA'].append(self.parse_hexa_fixed(l))
            elif l.startswith('CQUAD8'):
                elems['QUAD8'].append(self.parse_quad8_fixed(l))
            else: pass

        new_elems = {}
        if len(elems['TETRA']) > 0:
            TETRA = np.vstack([np.array((int(g[3]), int(g[4]), int(g[5]), int(g[6]),))
                           for g in elems['TETRA']])
            TETRA_ATTR = np.array([int(g[2]) for g in elems['TETRA']]) #PSOLID ID
            new_elems['TETRA'] = TETRA-1
            new_elems['TETRA_ATTR'] = TETRA_ATTR
            
        if len(elems['TRIA6']) > 0:                          
            TRIA6 = np.vstack([np.array((int(g[3]), int(g[4]), int(g[5]), ))
                           for g in elems['TRIA6']])
            TRIA6_ATTR = np.array([int(g[2]) for g in elems['TRIA6']])  #PSHELL ID
            new_elems['TRIA6'] = TRIA6-1
            new_elems['TRIA6_ATTR'] = TRIA6_ATTR
            
        if len(elems['HEXA']) > 0:                                  
            HEXA = np.vstack([np.array((int(g[3]), int(g[4]), int(g[5]), int(g[6]),
                           int(g[7]), int(g[8]), int(g[9]), int(g[10]),))
                           for g in elems['HEXA']])
            HEXA_ATTR = np.array([int(g[2]) for g in elems['HEXA']]) #PSOLID ID
            new_elems['HEXA'] = HEXA-1
            new_elems['HEXA_ATTR'] = HEXA_ATTR

        if len(elems['QUAD8']) > 0:
            QUAD8 = np.vstack([np.array((int(g[3]), int(g[4]), int(g[5]), int(g[6])))
                           for g in elems['QUAD8']])
            QUAD8_ATTR = np.array([int(g[2]) for g in elems['QUAD8']])  #PSHELL ID
            new_elems['QUAD8'] = QUAD8-1
            new_elems['QUAD8_ATTR'] = QUAD8_ATTR

        elems =  new_elems
        
        props = {'PSOLID':[],
                 'PSHELL':[]}
        while True:            
            l = self._read_line(fid)
            if l.startswith('ENDDATA'):break
            if l.startswith('PSOLID'):
                props['PSOLID'].append(self.parse_psolid_fixed(l))
            elif l.startswith('PSHELL'):
                props['PSHELL'].append(self.parse_pshell_fixed(l))
            else: pass
                          
        PSHELL = np.array([int(g[1]) for g in props['PSHELL']])  #PSHELL
        PSOLID = np.array([int(g[1]) for g in props['PSOLID']])  #PSOLID

        props = {'PSOLID': PSOLID, 
                 'PSHELL': PSHELL }


        dataset = {'PROPS':props,
                   'ELEMS' :elems,
                   'GRIDS':grids}
        fid.close()
        self.dataset = dataset

    def plot_tet(self, idx, **kwargs):
        from ifigure.interactive import solid
        
        grids = self.dataset['GRIDS']
        tet = self.dataset['ELEMS']['TETRA']
        i = tet[idx]
        pts = [grids[[i[0], i[1], i[2]]],
               grids[[i[1], i[2], i[3]]],
               grids[[i[2], i[3], i[0]]],
               grids[[i[3], i[0], i[1]]]]
        pts = np.rollaxis(np.dstack(pts), 2, 0)
        solid(pts, **kwargs)
        
    def _read_line(self, fid):
        cl = ''
        while True:
            line =  fid.readline()
            l = line.rstrip("\r\n")
            if l.startswith('+CONT'): l = cl + l[8:]
            if l.strip().endswith('+CONT'):
                cl = l.strip()[:-5]
                continue
            break
        return l
                
    def parse_grid_fixed(self, l):
        d=8
        cards= [l[d*i:d*(i+1)].strip() for i in range(8)]
        return cards
    def parse_tetra_fixed(self, l):
        d=8
        cards= [l[d*i:d*(i+1)].strip() for i in range(13)]
        return cards
    def parse_tria6_fixed(self, l):
        d=8
        cards= [l[d*i:d*(i+1)].strip() for i in range(9)]
        return cards
    def parse_hexa_fixed(self, l):
        d=8
        cards= [l[d*i:d*(i+1)].strip() for i in range(23)]
        return cards
    def parse_quad8_fixed(self, l):
        d=8
        cards= [l[d*i:d*(i+1)].strip() for i in range(11)]
        return cards        
    def parse_pshell_fixed(self, l):
        d=8
        cards= [l[d*i:d*(i+1)].strip() for i in range(2)]
        return cards
    def parse_psolid_fixed(self, l):
        d=8
        cards= [l[d*i:d*(i+1)].strip() for i in range(2)]
        return cards

def write_nas2mfem(filename,  reader, exclude_bdr = None):

        geom_type = {'TETRA': 4,
                     'TRIA6': 2,
                     'HEXA':  5,
                     'QUAD8': 3,}                     
        '''                     
        SEGMENT = 1
        TRIANGLE = 2
        SQUARE = 3
        TETRAHEDRON = 4
        CUBE = 5
        '''

        if exclude_bdr is None: exclude_bdr = [] 
        if reader.dataset is None:
            reader.load()
        data = reader.dataset
        fid = open(filename, 'w')

        grid = data['GRIDS']
        
        el_3d = ['TETRA','HEXA']
        el_2d = ['TRIA6', 'QUAD8']
        elems = data['ELEMS']
        unique_grids = list(np.unique(np.hstack([elems[name].flatten() for name in el_3d+el_2d if name in elems])))
        nvtc = len(unique_grids)
        ndim = grid.shape[-1]
        nelem = 0
        nbdry = 0        
        for k in elems:
            if k in el_3d: nelem = nelem + len(elems[k+'_ATTR'])
            if k in el_2d:
                tmp = [x for x in elems[k+'_ATTR'] if not x in exclude_bdr]
                nbrdy = nbdry + len(tmp)
        
        
        fid.write('MFEM mesh v1.0\n')
        fid.write('\n')
        fid.write('dimension\n')
        fid.write(str(ndim) + '\n')
        fid.write('\n')
        fid.write('elements\n')
        fid.write(str(nelem) + '\n')
        
        for name in el_3d:
            if not name in elems: continue            
            vidx = elems[name]
            attr = elems[name+'_ATTR']
            gtyp = geom_type[name]
            for i in range(len(attr)):

                txt = [str(attr[i]), str(gtyp)]
                txt.extend([str(unique_grids.index(x)) for x in vidx[i]])
                fid.write(' '.join(txt)+ '\n')
        fid.write('\n')                
        fid.write('boundary\n')
        fid.write(str(nbrdy) + '\n')
        for name in el_2d:
            if not name in elems: continue
            vidx = elems[name]
            attr = elems[name+'_ATTR']
            gtyp = geom_type[name]
            for i in range(len(attr)):
                if attr[i] in exclude_bdr: continue
                txt = [str(attr[i]), str(gtyp)]
                txt.extend([str(unique_grids.index(x)) for x in vidx[i]])
                fid.write(' '.join(txt) + '\n')
        fid.write('\n')                
        fid.write('vertices\n')
        fid.write(str(nvtc) + '\n')                        
        fid.write(str(ndim) + '\n')
        for i in range(nvtc):
            txt = [str(x) for x in grid[unique_grids[i]]]
            fid.write(' '.join(txt) + '\n')                          
        fid.close()

        
    
