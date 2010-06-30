import os.path
from optparse import OptionParser

def main():
    global OPTIONS
    usage = "usage: %prog [options] files"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--category", dest="category",
                        type='string', default="Anonymous", 
                        help="If provided, catagory name")
                       
    (OPTIONS, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("incorrect number of arguments")

    print '<Category name="%s">' % OPTIONS.category
    for arg in args:
        path = os.path.abspath(arg)
        name = os.path.basename(arg)
        print '  <Resource name="%s" path="%s" serverURL="http://www.broad.mit.edu/webservices/igv"/>' % (name, path)
    print '</Category>'

if __name__ == "__main__":
    main()
