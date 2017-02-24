from argparse import ArgumentParser

class ArgParser(ArgumentParser):
    def __init__(self, *args, **kwargs):
        self.argument_list = []
        ArgumentParser.__init__(self, *args, **kwargs)


    def add_argument(self, *args, **kwargs):
        ArgumentParser.add_argument(self, *args, **kwargs)
        self.argument_list.append(args[-1])

    def print_options(self, args):
        print("Options used:")
        d = vars(args)

        print self.argument_list
        keys = ['_'.join(filter(None, x.split('-'))) 
                for x in self.argument_list]
        print keys
        print d
        for k in keys:
            if k in d: print("   --" +  k + "  " + str(d[k]))

