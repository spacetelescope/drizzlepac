def findrootname(filename):
    puncloc = [filename.find(char) for char in string.punctuation]
    val = sys.maxint
    for num in puncloc:
        if num !=-1 and num < val:
            val = num
    return filename[0:val]
