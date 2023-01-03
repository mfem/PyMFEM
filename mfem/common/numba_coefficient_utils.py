'''
 subroutines for numba_coefficeint
'''


def generate_caller_scalar(setting):
    '''
    generate a callder function on the fly

    ex)
    if setting is
        {"isdepcomplex": (True, False), "kinds": (1, 0),
                       "output": True, size: (10, 1)}

    def _caller(ptx, sdim, data):
        ptx = farray(ptx, (sdim,), np.float64)      # for position
        arr0r = farray(data[0], (10,), np.float64)
        arr0i = farray(data[1], (10,), np.float64)
        arr0 = arr0r+1j*arr0i

        arr1 = farray(data[0], (1,), np.float64)

        params = (arr0, arr1)
        return (inner_func(ptx, *params))

    here inner_func is a function user provided.

    '''
    if setting['td']:
        text = ['def _caller(ptx, sdim, t, data):']
    else:
        text = ['def _caller(ptx, sdim, data):']

    text.append("    ptx = farray(ptx, (sdim,), np.float64)")
    count = 0
    params_line = '    params = ('

    for s, kind, size in zip(setting['isdepcomplex'], setting['kinds'], setting["sizes"]):
        if not isinstance(size, tuple):
            size = (size, )

        if s:
            t1 = '    arrr' + \
                str(count) + ' = farray(data[' + \
                str(count) + "], "+str(size) + ", np.float64)"
            t2 = '    arri' + \
                str(count) + ' = farray(data[' + \
                str(count+1) + "], "+str(size) + ", np.float64)"
            t3 = '    arr'+str(count) + ' = arrr' + \
                str(count) + "+1j*arri" + str(count)

            #if len(size) == 1 and size[0] == 1:
            if kind == 0:
                t1 += '[0]'
                t2 += '[0]'

            text.extend((t1, t2, t3))
            params_line += 'arr'+str(count)+','
            count = count + 2
        else:
            t = '    arr' + \
                str(count) + ' = farray(data[' + \
                str(count) + "], "+str(size) + ", np.float64)"

            #if len(size) == 1 and size[0] == 1:
            if kind == 0:            
                t += '[0]'

            text.append(t)

            params_line += 'arr'+str(count)+','
            count = count + 1

    params_line += ')'

    text.append(params_line)
    if setting["td"]:
        text.append("    return (inner_func(ptx, t, *params))")
    else:
        text.append("    return (inner_func(ptx, *params))")
    return '\n'.join(text)


def generate_caller_array_oldstyle(setting):
    '''
    generate a callder function on the fly
    ex)
    if setting is
        {"isdepcomplex": (True, False), "kinds": (1, 0),
                       "output": True, size: ((3, 3), 1), outsize: (2, 2) }
    def _caller(ptx, sdim, data, out_):
        ptx = farray(ptx, (sdim,), np.float64)      # for position
        arr0r = farray(data[0], (3, 3), np.float64)
        arr0i = farray(data[1], (3, 3), np.float64)
        arr0 = arr0r+1j*arr0i
        arr1 = farray(data[0], (1,), np.float64)
        out = farray(out_, (2, 2), np.complex128)
        params = (arr0, arr1, out,)
        return (inner_func(ptx, *params))
    here inner_func is a function user provided.
    '''
    if setting['td']:
        text = ['def _caller(ptx, sdim, t, data, out_):']
    else:
        text = ['def _caller(ptx, sdim, data, out_):']
    text.append("    ptx = farray(ptx, (sdim,), np.float64)")
    count = 0
    params_line = '    params = ('

    for s, kind, size in zip(setting['isdepcomplex'], setting['kinds'], setting["sizes"]):
        if not isinstance(size, tuple):
            size = (size, )

        if s:
            t1 = '    arrr' + \
                str(count) + ' = farray(data[' + \
                str(count) + "], "+str(size) + ", np.float64)"
            t2 = '    arri' + \
                str(count) + ' = farray(data[' + \
                str(count+1) + "], "+str(size) + ", np.float64)"

            #if len(size) == 1 and size[0] == 1:
            if kind == 0:                
                t1 += '[0]'
                t2 += '[0]'

            t3 = '    arr'+str(count) + ' = arrr' + \
                str(count) + "+1j*arri" + str(count)

            text.extend((t1, t2, t3))
            params_line += 'arr'+str(count)+','
            count = count + 2
        else:
            t = '    arr' + \
                str(count) + ' = farray(data[' + \
                str(count) + "], "+str(size) + ", np.float64)"

            #if len(size) == 1 and size[0] == 1:
            if kind == 0:                                        
                t += '[0]'

            text.append(t)

            params_line += 'arr'+str(count)+','
            count = count + 1

    outsize = setting["outsize"]
    if setting["output"]:
        t = '    out = farray(out_,' + str(outsize) + ", np.complex128)"
    else:
        t = '    out = farray(out_,' + str(outsize) + ", np.float64)"
    text.append(t)

    params_line += 'out, )'

    text.append(params_line)
    if setting["td"]:
        text.append("    return (inner_func(ptx, t, *params))")
    else:
        text.append("    return (inner_func(ptx, *params))")
    return '\n'.join(text)


def generate_caller_array(setting):
    '''
    generate a following callder function from setting

    ex)
    if setting is
        {"isdepcomplex": (True, False), "kinds": (1, 0),
                       "output": True, size: ((3, 3), 1), outsize: (2, 2) }

    def _caller(ptx, sdim, data, out_):
        ptx = farray(ptx, (sdim,), np.float64)      # for position
        arr0r = farray(data[0], (3, 3), np.float64)
        arr0i = farray(data[1], (3, 3), np.float64)
        arr0 = arr0r+1j*arr0i

        arr1 = farray(data[0], (1,), np.float64)

        out = farray(out_, (2, 2), np.complex128)

        params = (arr0, arr1, )

        ret = inner_func(ptx, *params)
        for i0 in range(2):
           for i1 in range(2):
              ret[i0,i1] = out[i0, i1]

    here inner_func is a function user provided.

    '''
    if setting['td']:
        text = ['def _caller(ptx, sdim, t, data, out_):']
    else:
        text = ['def _caller(ptx, sdim, data, out_):']
    text.append("    ptx = farray(ptx, (sdim,), np.float64)")
    count = 0
    params_line = '    params = ('

    for s, kind, size in zip(setting['isdepcomplex'], setting['kinds'], setting["sizes"]):
        if not isinstance(size, tuple):
            size = (size, )

        if s:
            t1 = '    arrr' + \
                str(count) + ' = farray(data[' + \
                str(count) + "], "+str(size) + ", np.float64)"
            t2 = '    arri' + \
                str(count) + ' = farray(data[' + \
                str(count+1) + "], "+str(size) + ", np.float64)"

            #if len(size) == 1 and size[0] == 1:
            if kind == 0:                                
                t1 += '[0]'
                t2 += '[0]'

            t3 = '    arr'+str(count) + ' = arrr' + \
                str(count) + "+1j*arri" + str(count)

            text.extend((t1, t2, t3))
            params_line += 'arr'+str(count)+','
            count = count + 2
        else:
            t = '    arr' + \
                str(count) + ' = farray(data[' + \
                str(count) + "], "+str(size) + ", np.float64)"

            #if len(size) == 1 and size[0] == 1:
            if kind == 0:                                                
                t += '[0]'

            text.append(t)

            params_line += 'arr'+str(count)+','
            count = count + 1

    outsize = setting["outsize"]
    if setting["output"]:
        t = '    out = farray(out_,' + str(outsize) + ", np.complex128)"
    else:
        t = '    out = farray(out_,' + str(outsize) + ", np.float64)"
    text.append(t)
    '''
    params_line += 'out, )'
    '''
    params_line += ')'
    text.append(params_line)
    if setting["td"]:
        text.append("    ret = inner_func(ptx, t, *params)")
    else:
        text.append("    ret = inner_func(ptx, *params)")

    idx_text = ""
    for k, s in enumerate(setting["outsize"]):
        text.append("    " + " "*k + "for i" + str(k) +
                    " in range(" + str(s) + "):")
        idx_text = idx_text + "i"+str(k)+","
    text.append("     " + " "*len(setting["outsize"]) +
                "out["+idx_text + "]=ret[" + idx_text + "]")

    return '\n'.join(text)


def generate_signature_scalar(setting):
    '''
    generate a signature to numba-compile a user scalar function

    ex)
    when user function is
        func(ptx, complex_array, float_scalar)

    setting is
        {"isdepcomplex": (2, 1), "kinds": (1, 0), "output": 2}

    output is
         types.complex128(types.double[:], types.complex128[:], types.double,)

    user function is

    '''
    sig = ''
    if setting['output']:
        sig += 'types.complex128(types.double[:], '
    else:
        sig += 'types.float64(types.double[:], '

    if setting['td']:
        sig += 'types.double, '

    for s, kind, in zip(setting['isdepcomplex'], setting['kinds'],):
        if s:
            if kind == 0:
                sig += 'types.complex128,'
            elif kind == 1:
                sig += 'types.complex128[:], '
            else:
                sig += 'types.complex128[:, :], '
        else:
            if kind == 0:
                sig += 'types.double,'
            elif kind == 1:
                sig += 'types.double[:], '
            else:
                sig += 'types.double[:, :], '

    sig = sig + ")"
    return sig


def generate_signature_array_oldstyle(setting):
    '''
    generate a signature to numba-compile a user scalar function

    ex)
    when user function is
        func(ptx, complex_array, float_scalar, complex_output_array_)

    setting is
        {"isdepcomplex": (2, 1), "kinds": (1, 0), "output": 2, "outkind": 2}

    output is
         types.void(types.double[:], types.complex128[::],
                    types.double, types.complex128[:, :])

    user function is

    '''
    sig = 'types.void(types.double[:], '
    if setting['td']:
        sig += 'types.double, '

    for s, kind, in zip(setting['isdepcomplex'], setting['kinds'],):
        if s:
            if kind == 0:
                sig += 'types.complex128,'
            elif kind == 1:
                sig += 'types.complex128[:], '
            else:
                sig += 'types.complex128[:,:], '
        else:
            if kind == 0:
                sig += 'types.double,'
            elif kind == 1:
                sig += 'types.double[:], '
            else:
                sig += 'types.double[:, :], '

    if setting['output']:
        if setting['outkind'] == 1:
            sig += 'types.complex128[:], '
        else:
            sig += 'types.complex128[:, :], '
    else:
        if setting['outkind'] == 1:
            sig += 'types.double[:], '
        else:
            sig += 'types.double[:, :], '
    sig = sig + ")"
    return sig


def generate_signature_array(setting):
    '''
    generate a signature to numba-compile a user scalar function

    ex)
    when user function is
        func(ptx, complex_array, float_scalar)

    setting is
        {"isdepcomplex": (2, 1), "kinds": (1, 0), "output": 2}

    output is
         types.complex128[:, :](types.double[:], types.complex128[:], types.double,)

    user function is

    '''
    sig = ''
    if setting['output']:
        if setting['outkind'] == 1:
            sig += 'types.complex128[:](types.double[:], '
        else:
            sig += 'types.complex128[:,:](types.double[:], '
    else:
        if setting['outkind'] == 1:
            sig += 'types.float64[:](types.double[:], '
        else:
            sig += 'types.float64[:,:](types.double[:], '
    if setting['td']:
        sig += 'types.double, '

    for s, kind, in zip(setting['isdepcomplex'], setting['kinds'],):
        if s:
            if kind == 0:
                sig += 'types.complex128,'
            elif kind == 1:
                sig += 'types.complex128[:], '
            else:
                sig += 'types.complex128[:, :], '
        else:
            if kind == 0:
                sig += 'types.double,'
            elif kind == 1:
                sig += 'types.double[:], '
            else:
                sig += 'types.double[:, :], '

    sig = sig + ")"
    return sig


def _process_dependencies(dependencies, setting):
    from mfem import mfem_mode
    if mfem_mode == 'serial':
        from mfem.ser import (Coefficient,
                              VectorCoefficient,
                              MatrixCoefficient,
                              IsNumbaCoefficient)
    else:
        from mfem.par import (Coefficient,
                              VectorCoefficient,
                              MatrixCoefficient,
                              IsNumbaCoefficient)

    iscomplex = []
    sizes = []
    kinds = []
    s_coeffs = []
    v_coeffs = []
    m_coeffs = []
    ns_coeffs = []
    nv_coeffs = []
    nm_coeffs = []
    for x in dependencies:
        if isinstance(x, tuple) and x[1] is not None:
            iscomplex.append(1)
            xx = x[0]
            if x[0] is None or x[1] is None:
                assert False, "dependency has to have both real imaginary parts defined"
            if isinstance(xx, VectorCoefficient):
                assert x[0].GetVDim() == x[1].GetVDim(
                ), "real and imaginary has to have the same vdim"
            if isinstance(xx, MatrixCoefficient):
                h1, w1 = x[0].GetHeight(), x[0].GetWidth()
                h2, w2 = x[1].GetHeight(), x[1].GetWidth()
                assert h1 == w1, "matrix must be square"
                assert h1 == h2, "real and imaginary has to have the same vdim"
                assert w1 == w2, "real and imaginary has to have the same vdim"
        else:
            if isinstance(x, tuple):
                # treat complex without imaginary part as real
                x = x[0]
            xx = x
            if IsNumbaCoefficient(xx):
                if xx.IsOutComplex():
                    iscomplex.append(2)
                else:
                    iscomplex.append(0)
            else:
                iscomplex.append(0)

            if isinstance(xx, MatrixCoefficient):
                assert xx.GetHeight() == xx.GetWidth(), "matrix must be square"

        if isinstance(xx, Coefficient):
            kinds.append(0)
            sizes.append(1)
            s_coeffs.append(xx)
            if iscomplex[-1] > 0:
                if IsNumbaCoefficient(xx) and not isinstance(x, tuple):
                    ns_coeffs.append(xx)
                else:
                    s_coeffs.append(x[1])

        elif isinstance(xx, VectorCoefficient):
            kinds.append(1)
            sizes.append(xx.GetVDim())
            v_coeffs.append(xx)
            if iscomplex[-1] > 0:
                if IsNumbaCoefficient(xx) and not isinstance(x, tuple):
                    nv_coeffs.append(xx)
                else:
                    v_coeffs.append(x[1])

        elif isinstance(xx, MatrixCoefficient):
            kinds.append(2)
            sizes.append((xx.GetHeight(), xx.GetWidth()))
            m_coeffs.append(xx)
            if iscomplex[-1] > 0:
                if IsNumbaCoefficient(xx) and not isinstance(x, tuple):
                    nm_coeffs.append(xx)
                else:
                    m_coeffs.append(x[1])

        else:
            assert False, "unknown coefficient type" + str(type(xx))

    setting['s_coeffs'] = s_coeffs
    setting['v_coeffs'] = v_coeffs
    setting['m_coeffs'] = m_coeffs
    setting['ns_coeffs'] = ns_coeffs
    setting['nv_coeffs'] = nv_coeffs
    setting['nm_coeffs'] = nm_coeffs

    return iscomplex, sizes, kinds


def get_setting(outsize, iscomplex=False, dependencies=None, td=False):
    setting = {}
    if iscomplex:
        setting['output'] = True
    else:
        setting['output'] = False

    isdepcomplex, sizes, kinds = _process_dependencies(dependencies, setting)

    setting['isdepcomplex'] = isdepcomplex
    setting['kinds'] = kinds
    setting['outsize'] = outsize
    if isinstance(outsize, tuple):
        setting['outkind'] = len(outsize)
    else:
        setting['outkind'] = 0
    setting['sizes'] = sizes
    setting['td'] = td

    return setting
