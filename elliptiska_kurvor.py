from math import sqrt, ceil
from typing import Optional

class Punkt:
    def __init__(self, x: int, y: int, kurva: ElliptiskKurva, är_oändlig=False):
        self.x = x
        self.y = y
        self.är_oändlig = är_oändlig
        self.kurva = kurva

    @classmethod
    def oändligheten(cls, kurva):
        return cls(-1, -1, kurva, är_oändlig=True)

    def __eq__(self, other):
        if not (type(other) is Punkt):
            return False
        if self.kurva != other.kurva:
            return False
        return (self.x == other.x) and (self.y == other.y) and (self.är_oändlig == other.är_oändlig)

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        # https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication

        if type(other) != Punkt:
            raise TypeError("Kan endast addera en punkt med en annan")
        if self.kurva != other.kurva:
            raise TypeError("Punkter som adderas måste tillhöra samma elliptiska kurva")

        # Addition med identitet
        if self.är_oändlig:
            return other
        if other.är_oändlig:
            return self

        # Vertikala linjer
        if self.x == other.x and self.y != other.y:
            return self.kurva.oändligheten

        if self == other:
            if self.y == 0:
                # Försök av att dividera med 0
                return self.kurva.oändligheten
            k = self.division(3 * self.x**2 + self.kurva.a, 2 * self.y)
        else:
            k = self.division(other.y - self.y, other.x - self.x)

        x3 = (k**2 - self.x - other.x) % self.kurva.p
        y3 = (k * (self.x - x3) - self.y) % self.kurva.p

        return self.kurva.vid(x3, y3)

    def __neg__(self):
        if self.är_oändlig: return self
        inv = self.kurva.vid(self.x, (-self.y) % self.kurva.p)
        assert inv is not None
        return inv

    def invertera(self): return self.__neg__()

    def modulär_invers(self, k: int):
        return pow(k, self.kurva.p - 2, self.kurva.p)

    def division(self, a: int, b: int):
        return (a * self.modulär_invers(b)) % self.kurva.p

    def __mul__(self, other: int):
        # TODO: Implementera negativ multiplikation
        # Double-and-add-algoritmen

        if type(other) != int:
            raise TypeError("Kan endast multiplicera en punkt med ett heltal")

        if other < 0:
            return self.invertera() * other

        mask = other
        tvåpotens = self
        summa = self.kurva.oändligheten
        while mask:
            if mask & 1 == 1:
                assert summa is not None
                summa += tvåpotens
            assert tvåpotens is not None
            tvåpotens = tvåpotens + tvåpotens
            mask >>= 1

        return summa

    def __repr__(self):
        if self.är_oändlig:
            return "Punkt(Oändligheten)"
        return f"Punkt({str(self.x)}, {str(self.y)})"


class ElliptiskKurva:
    def __init__(self, a: int, b: int, p: int):
        self.a = a
        self.b = b
        self.p = p
        self.diskriminant = self._beräkna_diskriminant() % p

        self.punkter: list[Punkt] = []
        self.oändligheten: Optional[Punkt] = None

        if not self._är_primtal(p):
            raise ValueError(f"p={str(p)} är inte ett primtal")
        if p == 2 or p == 3:
            raise ValueError(f"Karaktäristiken (char) K = {str(p)}. Kurvor med char = 2 | 3 kan inte användas")
        if self.diskriminant == 0:
            raise ValueError(f"Diskriminanten {str(self._beräkna_diskriminant())} ≡ 0 (mod {str(p)})")

        self._beräkna_punkter()

    def vid(self, x, y, fel=True):
        for punkt in self.punkter:
            if punkt.x == x and punkt.y == y:
                return punkt

        if fel:
            raise ValueError(f"Ingen punkt vid x: {str(x)}; y: {str(y)}")
        return None

    def __getitem__(self, item):
        if type(item) == tuple:
            x, y = item
            return self.vid(x, y)
        raise IndexError("Kan endast indexera en tupel med (x: int, y: int)")

    @staticmethod
    def _är_primtal(p: int) -> bool:
        for i in range(2, ceil(sqrt(p)) + 1):
            if p % i == 0:
                return False
        return True

    def _beräkna_diskriminant(self) -> int:
        return -16 * (4 * (self.a ** 3) + (self.b ** 2) * 27)

    def _beräkna_punkter(self):
        # TODO: använd coolare algoritm
        kvadrattabell: dict[int, list[int]] = {}
        for i in range(self.p):
            _sq = (i ** 2) % self.p
            if not _sq in kvadrattabell.keys():
                kvadrattabell[_sq] = []
            kvadrattabell[_sq].append(i)

        for x in range(self.p):
            hl = (x ** 3 + x * self.a + self.b) % self.p
            if hl in kvadrattabell.keys():
                for y in kvadrattabell[hl]:
                    self.punkter.append(Punkt(x, y, self))

        self.oändligheten = Punkt.oändligheten(self)
        self.punkter.append(self.oändligheten)

    def __repr__(self):
        return f"Elliptisk kurva `E` i y^2 = x^3 + {str(self.a) if self.a != 1 else ""}x" + \
               f"{" + " + str(self.b) if self.b != 0 else ""} med {str(len(self.punkter))} punkter, projekterad över F{str(self.p)}"
