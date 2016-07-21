from ifigure.interactive import figure

def open_meshviewer():
    viewer = figure(); viewer.threed(True)
    viewer.xlabel('x')
    viewer.ylabel('y')
    viewer.zlabel('z')
    w = viewer.attach_to_main()    
    return w
