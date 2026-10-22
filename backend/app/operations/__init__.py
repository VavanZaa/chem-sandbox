"""Reactor and process-model integration boundary."""

from .capabilities import operation_capabilities
from .validation import validate_operation

__all__ = ["operation_capabilities", "validate_operation"]
