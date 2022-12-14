'''
 subroutines for numba_coefficeint
'''


def generate_caller_scaler(settings):
    '''
    generate a callder function on the fly

    ex)
    if setting is {"input": (2, 1), "output": 2}

    def _caller(ptx, data):
        arr0 = data[0]+1j*data[0+1]
        arr2 = data[2]
        params = (arr0,arr2,)
        return (inner_func(ptx, *params))

    here inner_func is a function user provided.

    '''
    text = ['def _caller(ptx, data):']
    count = 0

    params_line = '    params = ('
    for s in settings["input"]:
        if s == 2:
            t = '    arr' + \
                str(count) + ' = data[' + str(count) + \
                    ']+1j*data[' + str(count) + '+1]'
            params_line += 'arr'+str(count)+','
            count = count + 2
        else:
            t = '    arr'+str(count) + ' = data[' + str(count) + ']'
            params_line += 'arr'+str(count)+','
            count = count + 1
        text.append(t)
    params_line += ')'

    text.append(params_line)
    text.append("    return (inner_func(ptx, *params))")
    return '\n'.join(text)


def generate_signature_scalar(setting):
    '''
    generate a signature to numba-compile a user scalar function

    ex)
    if setting is {"input": (2, 1), "output": 2}

    output : types.complex128(CPointer(types.double),
                              types.complex128,types.double,)

    '''

    sig = ''
    if setting['output'] == 1:
        sig += 'types.float64(CPointer(types.double, '
    else:
        sig += 'types.complex128(CPointer(types.double), '

    for s in setting['input']:
        if s == 1:
            sig += 'types.double,'
        else:
            sig += 'types.complex128,'

    sig = sig + ")"
    return sig


def generate_signature_array(setting):
    '''
    generate a signature to numba-compile a user vector/matrix function

    ex)
    if setting is {"input": (2, 1), "output": 2}

    output : types.void(CPointer(types.double), types.complex128, types.double,
                        CPointer(types.complex128))

    '''
    sig = ''
    sig += 'types.void(CPointer(types.double, '
    for s in setting['input']:
        if s == 1:
            sig += 'types.double,'
        else:
            sig += 'types.complex128,'

    if setting['output'] == 1:
        sig += 'CPointer(types.double), '
    else:
        sig += 'CPointer(types.complex128), '

    sig = sig + ")"
    return sig


def _process_dependencies(dependencies):
    from mfem import mfem_mode
    if mfem_mode == 'serial'
        import mfem.ser import (Coefficient
        import
                                VectorCoefficient,
                                MatrixCoefficient)
    else:
        import mfem.par import (Coefficient
        import
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
            if isinstacne(xx, VectorCoefficient):
                assert x[0].GetVdim() == x[1].GetVdim() "real and imaginary has to have the same vdim"
            if isinstacne(xx, MatrixCoefficient):
                assert x[0].Height() == x[0].Width() "matrix must be square"
                assert x[0].Height() == x[1].Height() "real and imaginary has to have the same vdim"
                assert x[0].Width() == x[1].Width() "real and imaginary has to have the same vdim"
        else:
            iscomplex.append(False)
            xx = x
            if isinstacne(xx, MatrixCoefficient):
                assert xx.Height() == xx.Width() "matrix must be square"

        if isinstacne(xx], Coefficient):
                kinds.append(0)
                sizes.append(1)
                s_coeffs.append(xx)
                if iscomplex[-1]:
                    s_coeffs.append(x[1])
        elif isinstacne(x, VectorCoefficient):
                kinds.append(1)
                sizes.append(xx.GetVdim())
                v_coeffs.append(xx)
                if iscomplex[-1]:
                    v_coeffs.append(x[1])
        elif isinstacne(xx, MatrixCoefficient):
                kinds.append(2)
                sizes.append(xx.GetHeight()*xx.GetWidth())
                m_coeffs.append(xx)
                if iscomplex[-1]:
                    m_coeffs.append(x[1])
        else:
            assert False, "unknown coefficient type" + str(type(xx))

    return iscomplex, sizes, kinds, s_coeffs, v_coeffs, m_coeffs

def get_setting(iscomplex=False, dependencies=None):
    setting = {}
    if iscomplex:
        setting['output'] = 2
    else:
        setting['output'] = 1

    iscomplex, sizes, kinds, s_coeffs, v_coeffs, m_coeffs = _process_dependencies(dependencies)

    setting['iscomplex'] = iscomplex
    setting['kinds'] = kinds
    setting['s_coeffs'] = s_coeffs
    setting['v_coeffs'] = v_coeffs
    setting['m_coeffs'] = m_coeffs
    return setting
