from __future__ import annotations

from typing import Any


class EvolutionOperator:
    def __init__(
        self,
        timeaxis: Any,
        dim: int | None = None,
        hamiltonian: Any = None,
        data: Any = None,
    ) -> None:

        self.timeaxis = timeaxis
        self.hamiltonian = hamiltonian

        pass

    def apply(opvec: Any) -> Any:
        """Apply the evolution operator to an operator or state vector

        Returns
        -------
        StateVectorEvolution or OperatorEvolution

        """
        pass

    def apply_left(opvec: Any) -> Any:
        """Apply the evolution operator to an operator or vector on the left

        Returns
        -------
        StateVectorEvolution or OperatorEvolution

        """
        pass
