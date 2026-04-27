from __future__ import annotations

from typing import Any

from .saveable import Saveable


class triangle(Saveable):
    """Class representing a symmetric matrix by a linear list"""

    def __init__(self, N: int = 0) -> None:
        self.N = N

    def get_empty_list(self) -> list[None]:
        return self.get_list(init=None)

    def get_list(self, init: Any = None) -> list[Any]:
        return [init] * (((self.N**2) - self.N) // 2 + self.N)

    def locate(
        self, i: int, j: int, transpose: bool = True, report_transpose: bool = False
    ) -> Any:

        if ((i >= self.N) or (j >= self.N)) or ((i < 0) or (j < 0)):
            raise Exception("Index out of range")

        trans = False
        if j > i:
            if transpose:
                trans = True
                pom = i
                i = j
                j = pom
            else:
                raise Exception("Index out of range (transposed is in range)")

        I = 0
        for m in range(1, i + 1):
            I -= m - 1
        I += j - i
        I += i * self.N

        if report_transpose:
            return trans, I
        return I

    def indices(self, I: int) -> None:
        pass


# FIXME: define triangle_list class which can store using triangle
