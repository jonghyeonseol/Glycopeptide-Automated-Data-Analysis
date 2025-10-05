"""
Report Builder
Provides fluent API for building reports
"""

from typing import Dict, Any, List
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from .base_report import BaseReport
from src.logger_config import get_logger

logger = get_logger(__name__)


class ReportBuilder:
    """
    Builder for composing reports from multiple sections

    Usage:
        report = (ReportBuilder("Analysis Summary")
            .add_section("Data Info", data_info_text)
            .add_section("Statistics", stats_text)
            .build()
            .save("summary.txt"))
    """

    def __init__(self, title: str = "Analysis Report"):
        """
        Initialize builder

        Args:
            title: Report title
        """
        self.title = title
        self.sections: List[Dict[str, Any]] = []

    def add_section(
        self,
        title: str,
        content: str,
        order: int = None
    ) -> 'ReportBuilder':
        """
        Add a section

        Args:
            title: Section title
            content: Section content
            order: Display order (auto-assigned if None)

        Returns:
            Self for chaining
        """
        if order is None:
            order = len(self.sections)

        self.sections.append({
            'title': title,
            'content': content,
            'order': order
        })

        return self

    def add_from_dict(self, data: Dict[str, str]) -> 'ReportBuilder':
        """
        Add sections from dictionary

        Args:
            data: Dictionary of {title: content}

        Returns:
            Self for chaining
        """
        for title, content in data.items():
            self.add_section(title, content)

        return self

    def build(self) -> 'ComposedReport':
        """
        Build the report

        Returns:
            ComposedReport instance
        """
        report = ComposedReport(self.title)

        for section in self.sections:
            report.add_section(
                title=section['title'],
                content=section['content'],
                order=section['order']
            )

        return report


class ComposedReport(BaseReport):
    """
    Report composed from builder

    This is a concrete implementation that simply
    uses the sections added via the builder.
    """

    def generate(self, data: Dict[str, Any] = None) -> str:
        """
        Generate report (sections already added)

        Args:
            data: Optional data (not used for composed reports)

        Returns:
            Rendered report
        """
        return self.render()
