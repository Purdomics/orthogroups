"""=================================================================================================


Michael Gribskov     05 March 2024
================================================================================================="""

class Dumb:

    def __init__(self):
        self.long = 'ausoihfjuhwpoeijpo;fisajdofikjlkgfdja[soifjasoiudjfoisdufpoiauf'
        self.short = ''

    def __next__(self):
        return self

    def __iter__(self):
        pos = 0
        buffer = ''
        if buffer:
            yield self
        while pos < len(self.long) - 10:
            # yield self.short
            self.short = self.long[pos:pos + 10]
            yield self
            pos += 10

        # raise StopIteration

        return

    def adv(self):
        return self.__iter__()
# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    a = Dumb()
    aiter = iter(a)
    for aa in aiter:
        # print(f'a:{a.short}\t{aa}\t{aa.short}')
        b = next(aiter)
        print(f'a:{a.short}\t{aa}')


    exit(0)
