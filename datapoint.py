class DataPoint:
    def __init__(self, t=0, k=1, cp=1, p=4, s=False):
        self.__k = k
        self.__t = t
        self.__cp = cp
        self.__p = p
        self.__static = s

    def set_t(self, temp):
        if not self.__static:
            self.__t = temp

    def get_p(self):
        return self.__p

    def get_cp(self):
        return self.__cp

    def get_t(self):
        return self.__t

    def get_k(self):
        return self.__k
    
    def __repr__(self):
        return str(self.__t)


