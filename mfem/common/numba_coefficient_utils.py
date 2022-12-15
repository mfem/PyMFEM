'''
 subroutines for numba_coefficeint
'''


def generate_caller_scalar(settings):
    '''
    generate a callder function on the fly

    ex)
    if setting is
        {"iscomplex": (True, False), "kinds": (1, 0),
                       "output": True, size: (10, 1)}

    def _caller(ptx, data):
        ptx = numba.carray(ptx, (sdim,), np.float64)      # for position
        arr0r = numba.carray(data[0], (10,), np.float64)
        arr0i = numba.carray(data[1], (10,), np.float64)
        arr0 = arr0r+1j*arr0i

        arr1 = numba.carray(data[0], (1,), np.float64)

        params = (arr0, arr1)
        return (inner_func(ptx, *params))

    here inner_func is a function user provided.

    '''
    text = ['def _caller(ptx, data):']
    text.append("    ptx = numba.carray(ptx, (sdim,), np.float64)")
    count = 0
    params_line = '    params = ('

    for s, kind, size in zip(setting['iscomplex'], setting['kinds'], setting("size")):
        if s:
            t1 = '    arrr' + \
                str(count) + ' = numba.carray[data[' + \
                str(count) + "], ("+str(size) + "), np.float64)"
            t2 = '    arri' + \
                str(count) + ' = numba.carray[data[' + \
                str(count+1) + "], ("+str(size) + "), np.float64)"
            t3 = '    arr'+str(count) + ' = arrr' + \
                str(count) + "+1j*arri" + str(count)

            text.extend((t1, t2, t3))
            params_line += 'arr'+str(count)+','
            count = count + 2
        else:
            t = '    arr' + \
                str(count) + ' = numba.carray[data[' + \
                str(count) + "], ("+str(size) + "), np.float64)"
            text.append(t)

            params_line += 'arr'+str(count)+','
            count = count + 1

    params_line += ')'

    text.append(params_line)
    text.append("    return (inner_func(ptx, *params))")
    return '\n'.join(text)


def generate_caller_array(settings):
    '''
    generate a callder function on the fly

    ex)
    if setting is
        {"iscomplex": (True, False), "kinds": (1, 0),
                       "output": True, size: ((3, 3), 1), outsize: (2, 2) }

    def _caller(ptx, data):
        ptx = numba.carray(ptx, (sdim,), np.float64)      # for position
        arr0r = numba.carray(data[0], (3, 3), np.float64)
        arr0i = numba.carray(data[1], (3, 3), np.float64)
        arr0 = arr0r+1j*arr0i

        arr1 = numba.carray(data[0], (1,), np.float64)

        outr = numba.carray(data[0], (2, 2), np.float64)
        outi = numba.carray(data[1], (2, 2), np.float64)
        out = outr+1j*outi

        params = (arr0, arr1, out)
        return (inner_func(ptx, *params))

    here inner_func is a function user provided.

    '''
    text = ['def _caller(ptx, data):']
    text.append("    ptx = numba.carray(ptx, (sdim,), np.float64)")
    count = 0
    params_line = '    params = ('

    for s, kind, size in zip(setting['iscomplex'], setting['kinds'], setting("size")):
        if s:
            if not isinstances(size, tuple):
                size = (size, )
            t1 = '    arrr' + \
                str(count) + ' = numba.carray[data[' + \
                str(count) + "], "+str(size) + ", np.float64)"
            t2 = '    arri' + \
                str(count) + ' = numba.carray[data[' + \
                str(count+1) + "], "+str(size) + ", np.float64)"
            t3 = '    arr'+str(count) + ' = arrr' + \
                str(count) + "+1j*arri" + str(count)

            text.extend((t1, t2, t3))
            params_line += 'arr'+str(count)+','
            count = count + 2
        else:
            t = '    arr' + \
                str(count) + ' = numba.carray[data[' + \
                str(count) + "], ("+str(size) + "), np.float64)"
            text.append(t)

            params_line += 'arr'+str(count)+','
            count = count + 1

    if output:
        t1 = '    outr = numba.carray[data[' + \
            str(count) + "], "+str(outsize) + ", np.float64)"
        t2 = '    outi = numba.carray[data[' + \
            str(count+1) + "], "+str(outsize) + ", np.float64)"
        t3 = '    out = outr+1j*outi'
        text.extend((t1, t2, t3))
    else:
        t1 = '    out = numba.carray[data[' + \
            str(count) + "], "+str(outsize) + ", np.float64)"
        text.append(t)

    params_line += 'out)'

    text.append(params_line)
    text.append("    return (inner_func(ptx, *params))")
    return '\n'.join(text)


def generate_signature_scalar(setting):
    '''
    generate a signature to numba-compile a user scalar function

    ex)
    when user function is
        func(ptx, complex_array, float_scalar)

    setting is
        {"iscomplex": (2, 1), "kinds": (1, 0), "output": 2}

    output is
         types.complex128(types.double[:], types.complex128[:], types.double,)

    user function is

    '''
    sig = ''
    if setting['output'] == 1:
        sig += 'types.float64(types.double[:], '
    else:
        sig += 'types.complex128(types.double[:], '

    for s, kind, in zip(setting['iscomplex'], setting['kinds'],):
        if s:
            if kind == 0:
                sig += 'types.complex128,'
            else:
                sig += 'types.complex128[:], '
        else:
            if kind == 0:
                sig += 'types.double,'
            else:
                sig += 'types.double[:], '

    sig = sig + ")"
    return sig


def generate_signature_array(setting):
    '''
    generate a signature to numba-compile a user scalar function

    ex)
    when user function is
        func(ptx, complex_array, float_scalar, complex_output_array_)

    setting is
        {"iscomplex": (2, 1), "kinds": (1, 0), "output": 2}

    output is
         types.void(types.double[:], types.complex128[:],
                    types.double, types.complex128[:])

    user function is

    '''
    sig = 'types.void(types.double[:], '
    for s, kind, in zip(setting['iscomplex'], setting['kinds'],):
        if s:
            if kind == 0:
                sig += 'types.complex128,'
            else:
                sig += 'types.complex128[:], '
        else:
            if kind == 0:
                sig += 'types.double,'
            else:
                sig += 'types.double[:], '

    if setting['output']:
        sig += 'types.complex128[:], '
    else:
        sig += 'types.double[:], '

    sig = sig + ")"
    return sig


def _process_dependencies(dependencies):
    from mfem import mfem_mode
    if mfem_mode == 'serial':
        from mfem.ser import (Coefficient,
                              VectorCoefficient,
                              MatrixCoefficient)
    else:
        from mfem.par import (Coefficient,
                              VectorCoefficient,
                              MatrixCoefficient)
    iscomplex = []
    sizes = []
    kinds = []
    s_coeffs = []
    v_coeffs = []
    m_coeffs = []
    for x in dependencies:
        if isinstance(x, tuple):
            iscomplex.append(True)
            xx = x[0]
            if x[0] is None or x[1] is None:
                assert False, "dependency has to have both real imaginary parts defined"
            if isinstance(xx, VectorCoefficient):
                assert x[0].GetVdim() == x[1].GetVdim(
                ), "real and imaginary has to have the same vdim"
            if isinstance(xx, MatrixCoefficient):
                h1, w1 = x[0].Height(), x[0].Width()
                h2, w2 = x[1].Height(), x[1].Width()
                assert h1 == w1, "matrix must be square"
                assert h1 == h2, "real and imaginary has to have the same vdim"
                assert w1 == w2, "real and imaginary has to have the same vdim"
        else:
            iscomplex.append(False)
            xx = x
            if isinstance(xx, MatrixCoefficient):
                assert xx.Height() == xx.Width(), "matrix must be square"

        if isinstance(xx, Coefficient):
            kinds.append(0)
            sizes.append(1)
            s_coeffs.append(xx)
            if iscomplex[-1]:
                s_coeffs.append(x[1])
        elif isinstance(x, VectorCoefficient):
            kinds.append(1)
            sizes.append(xx.GetVdim())
            v_coeffs.append(xx)
            if iscomplex[-1]:
                v_coeffs.append(x[1])
        elif isinstance(xx, MatrixCoefficient):
            kinds.append(2)
            sizes.append(xx.GetHeight()*xx.GetWidth())
            m_coeffs.append(xx)
            if iscomplex[-1]:
                m_coeffs.append(x[1])
        else:
            assert False, "unknown coefficient type" + str(type(xx))

    return iscomplex, sizes, kinds, s_coeffs, v_coeffs, m_coeffs


def get_setting(outsize, iscomplex=False, dependencies=None):
    setting = {}
    if iscomplex:
        setting['output'] = True
    else:
        setting['output'] = False

    iscomplex, sizes, kinds, s_coeffs, v_coeffs, m_coeffs = _process_dependencies(
        dependencies)

    setting['iscomplex'] = iscomplex
    setting['kinds'] = kinds
    setting['s_coeffs'] = s_coeffs
    setting['v_coeffs'] = v_coeffs
    setting['m_coeffs'] = m_coeffs
    setting['outsize'] = outsize

    return setting
