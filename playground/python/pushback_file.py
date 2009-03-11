class pushback_file(file):
    """Opens a file using the standard file interface adding the ability 
to pushback or unread some section of the file that was read from the file."""

    def __init__(self, fname, mode='r', bufsize=0):
        file.__init__(self, fname, mode, bufsize)
        self.pushed_back = []
        
    def next(self):
        if len(self.pushed_back):
            return self.pushed_back.pop()
        else:
            return file.next(self)
        
    def pushback(self, item):
        """Put some item (bytes, line, etc.) back on the \"front\" of the file"""
        self.pushed_back.append(item)
