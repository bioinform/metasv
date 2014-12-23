class Enum(set):
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

SV_Types = Enum(["DEL", "INS", "INV", "DUP", "DUP:TANDEM", "ITX", "CTX"])
