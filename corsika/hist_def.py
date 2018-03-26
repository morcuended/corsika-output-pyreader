class Histogram:
    pass


__m = Histogram()  # m will contain all module-level values
__m.dbname = None  # database name global in module

# later, in code:
if __m.dbname is None:
    __m.dbname = name