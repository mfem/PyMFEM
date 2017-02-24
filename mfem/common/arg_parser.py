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

        keys = ['_'.join(filter(None, x.split('-'))) 
                for x in self.argument_list]
        for k in keys:
            if k in d: print("   --" +  k + "  " + str(d[k]))

