"""
Workflow Definitions
Provides composable workflow steps for pipeline execution
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Callable, Optional
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class WorkflowStatus(Enum):
    """Workflow execution status"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


class WorkflowStep(ABC):
    """
    Abstract base class for workflow steps

    Each step performs a specific operation in the pipeline
    """

    def __init__(self, name: str, description: str = ""):
        self.name = name
        self.description = description
        self.status = WorkflowStatus.PENDING

    @abstractmethod
    def execute(self, state: Any, config: Dict[str, Any]) -> None:
        """
        Execute the workflow step

        Args:
            state: Pipeline state object
            config: Configuration dictionary
        """

    def on_start(self) -> None:
        """Hook called before step executes"""
        self.status = WorkflowStatus.RUNNING
        logger.info(f"\n[{self.name}] Starting...")
        if self.description:
            logger.info(f"  Description: {self.description}")

    def on_complete(self) -> None:
        """Hook called after step completes"""
        self.status = WorkflowStatus.COMPLETED
        logger.info(f"[{self.name}] Completed ✓")

    def on_error(self, error: Exception) -> None:
        """Hook called when step fails"""
        self.status = WorkflowStatus.FAILED
        logger.error(f"[{self.name}] Failed: {str(error)}")

    def run(self, state: Any, config: Dict[str, Any]) -> None:
        """
        Run the workflow step with hooks

        Args:
            state: Pipeline state
            config: Configuration
        """
        try:
            self.on_start()
            self.execute(state, config)
            self.on_complete()
        except Exception as e:
            self.on_error(e)
            raise


class Workflow:
    """
    Workflow composed of multiple steps

    Executes steps in sequence, managing state between them
    """

    def __init__(self, name: str, steps: Optional[list] = None):
        """
        Initialize workflow

        Args:
            name: Workflow name
            steps: List of WorkflowStep instances
        """
        self.name = name
        self.steps = steps or []
        self.status = WorkflowStatus.PENDING

    def add_step(self, step: WorkflowStep) -> 'Workflow':
        """
        Add a step to the workflow

        Args:
            step: WorkflowStep to add

        Returns:
            Self for chaining
        """
        self.steps.append(step)
        return self

    def execute(self, state: Any, config: Dict[str, Any]) -> None:
        """
        Execute all steps in the workflow

        Args:
            state: Pipeline state
            config: Configuration
        """
        logger.info(f"\n{'=' * 80}")
        logger.info(f"Executing Workflow: {self.name}")
        logger.info(f"{'=' * 80}")

        self.status = WorkflowStatus.RUNNING

        try:
            for i, step in enumerate(self.steps, 1):
                logger.info(f"\nStep {i}/{len(self.steps)}")
                step.run(state, config)

            self.status = WorkflowStatus.COMPLETED
            logger.info(f"\nWorkflow '{self.name}' completed successfully\n")

        except Exception as e:
            self.status = WorkflowStatus.FAILED
            logger.error(f"Workflow '{self.name}' failed: {str(e)}")
            raise


class FunctionalStep(WorkflowStep):
    """
    Workflow step that wraps a function

    Allows creating steps from functions without subclassing
    """

    def __init__(self, name: str, func: Callable, description: str = ""):
        super().__init__(name, description)
        self.func = func

    def execute(self, state: Any, config: Dict[str, Any]) -> None:
        """Execute the wrapped function"""
        self.func(state, config)
