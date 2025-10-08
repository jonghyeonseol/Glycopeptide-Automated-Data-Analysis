"""
Base Report Module
Provides abstract base class for report generation
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any
from pathlib import Path
from dataclasses import dataclass

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.logger_config import get_logger

logger = get_logger(__name__)


@dataclass
class ReportSection:
    """
    Report section container

    Attributes:
        title: Section title
        content: Section content (text)
        order: Display order (lower = earlier)
    """
    title: str
    content: str
    order: int = 0


class BaseReport(ABC):
    """
    Abstract base class for reports

    Provides:
    - Section management
    - Template rendering
    - File saving
    - Consistent formatting
    """

    def __init__(self, title: str = "Analysis Report"):
        """
        Initialize base report

        Args:
            title: Report title
        """
        self.title = title
        self.sections: List[ReportSection] = []

    def add_section(
        self,
        title: str,
        content: str,
        order: int = 0
    ) -> 'BaseReport':
        """
        Add a section to the report

        Args:
            title: Section title
            content: Section content
            order: Display order

        Returns:
            Self for chaining
        """
        section = ReportSection(title=title, content=content, order=order)
        self.sections.append(section)
        return self

    @abstractmethod
    def generate(self, data: Dict[str, Any]) -> str:
        """
        Generate report content

        Args:
            data: Data for report generation

        Returns:
            Formatted report content
        """

    def render(self) -> str:
        """
        Render all sections into final report

        Returns:
            Complete report text
        """
        # Sort sections by order
        sorted_sections = sorted(self.sections, key=lambda s: s.order)

        # Build report
        lines = []
        lines.append("=" * 80)
        lines.append(self.title)
        lines.append("=" * 80)
        lines.append("")

        for section in sorted_sections:
            if section.title:
                lines.append(f"\n{section.title}")
                lines.append("-" * len(section.title))
            lines.append(section.content)
            lines.append("")

        return "\n".join(lines)

    def save(self, output_path: Path) -> Path:
        """
        Save report to file

        Args:
            output_path: Path to save report

        Returns:
            Path to saved file
        """
        content = self.render()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            f.write(content)

        logger.info(f"Saved report to {output_path}")
        return output_path

    def clear(self) -> None:
        """Clear all sections"""
        self.sections.clear()
