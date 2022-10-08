import os
path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

mfem_mode = None
pymfem_debug = -1

def debug_print(message):
    if pymfem_debug < 0:   # debug < 0 
        return
    elif pymfem_debug == 0:   # debug = 0
        pass
    elif pymfem_debug > 0:   # debug = 1
        pass
    elif pymfem_debug > 1: # debug = 2
        import traceback as tb
        tb.print_stack()
    elif pymfem_debug > 2: # debug = 3
        import traceback as tb
        print(''.join(tb.format_list([tb.extract_stack(limit = 2)[0],])))

    print(message)

__version__ = '4.4.0.1'

