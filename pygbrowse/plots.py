"""Browser tracks and the GenomeBrowser that stacks them into one figure.

Every track draws a 0-based half-open window ``[ws, we)`` of one chromosome (UCSC names)
on a matplotlib Axes. Glyph sizes given in inches are converted to data units with the
figure width and row height that ``GenomeBrowser.visualize`` passes to ``plot``.
"""

import abc
import logging
from collections.abc import Callable, Iterable, Sequence

import intervaltree
import matplotlib
import matplotlib.axes
import matplotlib.cm
import matplotlib.colors
import matplotlib.figure
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy
import pandas
from scipy import ndimage

from . import _arcs, utilities
from .datasources import IntervalData

logger = logging.getLogger(__name__)

DEFAULT_YLABEL_PAD = 50

# Track style, formerly applied globally through seaborn.set_style("ticks").
DEFAULT_STYLE = {
    "figure.facecolor": "white",
    "axes.labelcolor": ".15",
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.color": ".15",
    "ytick.color": ".15",
    "axes.axisbelow": True,
    "grid.linestyle": "-",
    "text.color": ".15",
    "font.family": ["sans-serif"],
    "font.sans-serif": [
        "Arial",
        "DejaVu Sans",
        "Liberation Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ],
    "lines.solid_capstyle": "round",
    "patch.edgecolor": "w",
    "patch.force_edgecolor": True,
    "xtick.top": False,
    "ytick.right": False,
    "axes.grid": False,
    "axes.facecolor": "#EAEAF2",
    "axes.edgecolor": "w",
    "grid.color": ".8",
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "axes.spines.right": True,
    "axes.spines.top": True,
    "xtick.bottom": True,
    "ytick.left": True,
}


class _BrowserSubPlot(abc.ABC):
    """Base class for one browser track."""

    @abc.abstractmethod
    def plot(
        self,
        ax: matplotlib.axes.Axes,
        chrom: str,
        ws: int,
        we: int,
        fig_width: float,
        row_height: float,
    ) -> object:
        """Draw the track.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes of this track's row.
        chrom : str
            Chromosome, UCSC style.
        ws, we : int
            Window, 0-based half-open.
        fig_width : float
            Figure width in inches.
        row_height : float
            Row height in inches.
        """


class InteractionPlot(_BrowserSubPlot):
    def __init__(
        self,
        data: pandas.DataFrame,
        bin_size: int,
        arc_color: object = (0.7, 0.3, 0.6),
        direction: str = "down",
        baseline: float | None = None,
        vertical_scaling_factor: float = 1,
        thickness_column: str | None = None,
        show_bin_centers: bool = True,
        label: str = "Plac-seq",
    ):
        """Draw pairwise interactions (e.g. PLAC-seq, Hi-C loops) as arcs.

        Parameters
        ----------
        data : pandas.DataFrame
            One row per interaction, with columns ``chr1``, ``start1``, ``end1``,
            ``start2``, ``end2`` (0-based half-open anchors).
        bin_size : int
            Anchor bin size in bp; used to place the bin-center ticks.
        arc_color : matplotlib color
            Arc color.
        direction : str
            ``down`` hangs arcs from the top of the row, ``up`` raises them from the bottom.
        baseline : float, optional
            y position of the anchors. Defaults to 1 for ``down`` and 0 for ``up``.
        vertical_scaling_factor : float
            Multiplier on arc heights.
        thickness_column : str, optional
            Reserved; not used yet.
        show_bin_centers : bool
            Draw an x tick at every bin center.
        label : str
            y-axis label.
        """
        self.data = data
        self.bin_size = bin_size
        if baseline is None:
            baseline = 1 if direction == "down" else 0
        self.baseline = baseline
        self.vertical_scaling_factor = vertical_scaling_factor
        self.arc_color = arc_color
        self.label = label
        self.direction = direction
        self.show_bin_centers = show_bin_centers
        self.thickness_column = thickness_column

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        """Draw every interaction with at least one anchor midpoint inside the window."""
        on_chrom = self.data.loc[self.data["chr1"] == chrom]
        left_midpoints = (on_chrom["end1"] + on_chrom["start1"]) / 2
        right_midpoints = (on_chrom["end2"] + on_chrom["start2"]) / 2
        left_visible = (left_midpoints >= ws) & (left_midpoints <= we)
        right_visible = (right_midpoints >= ws) & (right_midpoints <= we)
        visible_ids = on_chrom.index[left_visible | right_visible]

        original_ylim = ax.get_ylim()
        ax.set_xlim(ws, we)
        for interaction_id in visible_ids:
            _arcs.draw_arc_interaction(
                ax,
                left_bin_center=left_midpoints.loc[interaction_id],
                right_bin_center=right_midpoints.loc[interaction_id],
                xlim=(ws, we),
                ylim=original_ylim,
                color=self.arc_color,
                baseline=self.baseline,
                vertical_scaling_factor=self.vertical_scaling_factor,
                direction=self.direction,
            )
        ax.set_xlim(ws, we)
        ax.set_ylim(original_ylim)

        if self.label:
            ax.set_ylabel(self.label)
        if self.show_bin_centers:
            self._draw_bin_center_ticks(ax, ws, we)
        ax.set_yticks([])

    def _draw_bin_center_ticks(self, ax: matplotlib.axes.Axes, ws: int, we: int) -> None:
        """Put an unlabeled x tick at the center of every bin in the window."""
        half_bin = self.bin_size / 2
        leftmost_tick = numpy.ceil((ws - half_bin) / self.bin_size) * self.bin_size + half_bin
        rightmost_tick = numpy.floor((we - half_bin) / self.bin_size + 1) * self.bin_size + half_bin
        ax.set_xticks(numpy.arange(leftmost_tick, rightmost_tick, self.bin_size))
        ax.set_xticklabels([])
        if self.direction == "down":
            ax.xaxis.set_ticks_position("top")


class BedPlot(_BrowserSubPlot):
    DEFAULT_PATCH_KWARGS = {"linewidth": 1, "edgecolor": "k"}

    def __init__(
        self,
        data: IntervalData,
        label: str = "",
        cmap: str = "RdBu_r",
        baseline: float = 0.5,
        color: object = "k",
        color_by: str = "",
        display_value: str = "",
        patch_height: float = 0.5,
        pad_fraction: float = 0.1,
        patch_kwargs: dict | None = None,
    ):
        """Draw intervals as boxes on one line of a row.

        Parameters
        ----------
        data : IntervalData
            Intervals to draw.
        label : str
            y tick label at ``baseline``. Several BedPlots can share one row.
        cmap : str
            Colormap for ``color_by``, centered on 0.
        baseline : float
            y position of the boxes' centers.
        color : matplotlib color
            Box color when ``color_by`` is not set.
        color_by : str
            Column whose value colors each box.
        display_value : str
            Column whose value is written on each box.
        patch_height : float
            Box height in data units.
        pad_fraction : float
            y padding fraction used when widening the y limits to fit the boxes.
        patch_kwargs : dict, optional
            Extra ``Rectangle`` keyword arguments, merged over ``DEFAULT_PATCH_KWARGS``.

        Raises
        ------
        ValueError
            If ``color_by`` or ``display_value`` is not a column of the interval data.
        """
        self.data = data
        self.color = color
        self.pad_fraction = pad_fraction

        columns = data.data.columns
        self.color_by = color_by
        if self.color_by:
            if self.color_by not in columns:
                raise ValueError(f"Color-by column {self.color_by!r} is not in the interval data")
            extent = numpy.abs(data.data[self.color_by]).max()
            self.color_mapper = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(vmin=-extent, vmax=extent), cmap=cmap
            )
        else:
            self.color_mapper = None

        self.display_value = display_value
        if self.display_value and self.display_value not in columns:
            raise ValueError(
                f"Display value column {self.display_value!r} is not in the interval data"
            )
        self.label = label
        self.patch_height = patch_height
        self.baseline = baseline
        self.patch_kwargs = {**self.DEFAULT_PATCH_KWARGS, **(patch_kwargs or {})}

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        """Draw every interval overlapping the window."""
        utilities.add_label(ax=ax, tick=self.baseline, tick_label=self.label, axis="y")

        for _, interval in self.data.query(chrom, ws, we).iterrows():
            self._draw_interval(ax, interval)

        for edge in (self.baseline + self.patch_height / 2, self.baseline - self.patch_height / 2):
            utilities.adjust_limits(
                ax=ax, new_position=edge, axis="y", padding_fraction=self.pad_fraction
            )

    def _draw_interval(self, ax: matplotlib.axes.Axes, interval: pandas.Series) -> None:
        """Draw one interval's box and, if requested, its value."""
        start_loc, end_loc = interval["chromStart"], interval["chromEnd"]
        if end_loc <= start_loc:
            raise ValueError(f"Interval end {end_loc} is not greater than start {start_loc}")
        if self.color_mapper is not None:
            interval_color = self.color_mapper.to_rgba(interval[self.color_by])
        else:
            interval_color = self.color

        ax.add_patch(
            matplotlib.patches.Rectangle(
                xy=(start_loc, self.baseline - self.patch_height / 2),
                width=end_loc - start_loc,
                height=self.patch_height,
                facecolor=interval_color,
                **self.patch_kwargs,
            )
        )
        if self.display_value:
            ax.text(
                x=(start_loc + end_loc) / 2,
                y=self.baseline,
                s=f"{interval[self.display_value]:>0.2}",
                ha="center",
            )


class WigPlot(_BrowserSubPlot):
    def __init__(
        self,
        data,
        label: str | None = None,
        color: object = None,
        solid: bool = True,
        alpha: float = 1.0,
        center_vector: bool = False,
        scale_vector_to_plot: bool = False,
        label_rotation: float = 0,
    ):
        """Draw a signal vector as a filled area or a line.

        Parameters
        ----------
        data : _VectorDataSource
            Signal source, e.g. ``BigWigDataSource``, ``TagDirectory`` or ``SparseVectors``.
            Smoothing and transforms are options of the data source.
        label : str, optional
            y-axis label and legend label.
        color : matplotlib color, optional
            Fill or line color.
        solid : bool
            Fill under the signal; otherwise draw a line.
        alpha : float
            Opacity.
        center_vector : bool
            Shift the signal so its mean sits at the middle of the row.
        scale_vector_to_plot : bool
            Scale the signal's range to the row's y range.
        label_rotation : float
            y-label rotation in degrees.
        """
        self.data = data
        self.color = color
        self.solid = solid
        self.alpha = alpha
        self.center = center_vector
        self.scale_vector_to_plot = scale_vector_to_plot
        self.label = label
        self.label_rotation = label_rotation

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        """Draw the signal over the window and autoscale the y axis."""
        ylim = ax.get_ylim()
        vert_span = ylim[1] - ylim[0]
        vert_center = vert_span / 2 + ylim[0]

        # astype returns a copy, so scaling below never writes into the data source.
        signal = self.data.query(chrom, ws, we).astype(float)
        if self.scale_vector_to_plot:
            signal = signal / (signal.max() - signal.min()) * vert_span
        if self.center:
            signal = signal - signal.mean() + vert_center
        signal = signal.loc[(signal.index >= ws) & (signal.index < we)]
        signal.name = self.label

        style = {"color": self.color, "alpha": self.alpha, "label": self.label}
        if self.solid:
            ax.fill_between(x=signal.index, y1=signal, **style)
        else:
            ax.plot(signal.index, signal, **style)
        ax.autoscale(enable=True, axis="y")

        if self.label:
            ax.set_ylabel(self.label, rotation=self.label_rotation, labelpad=DEFAULT_YLABEL_PAD)


class GeneModelPlot(_BrowserSubPlot):
    def __init__(
        self,
        data,
        label: str = "Genes",
        color: object = "k",
        feature_height: float = 0.12,
        chevron_height: float = 0.05,
        chevron_width: float = 0.04,
        chevron_spacing: float = 0.10,
        truncation_size: float = 0.10,
        utr_endcap_width: float = 0.04,
        gene_name_fontsize: float = 8,
        genes_include: Iterable[str] | None = None,
    ):
        """Draw gene models: backbone, strand chevrons, CDS boxes, UTRs and names.

        Overlapping genes are stacked on separate levels. Sizes are in inches.

        Parameters
        ----------
        data : _GeneModels
            Gene annotation source, e.g. ``Gff3Annotations``.
        label : str
            y-axis label.
        color : matplotlib color
            Color of all glyphs.
        feature_height : float
            Half-height of a CDS box (a UTR body is half as tall).
        chevron_height, chevron_width, chevron_spacing : float
            Strand chevron size and spacing.
        truncation_size : float
            Size of the arrowheads marking a gene cut off by the window edge.
        utr_endcap_width : float
            Width of the full-height cap on the CDS side of a UTR.
        gene_name_fontsize : float
            Gene name font size in points.
        genes_include : iterable of str, optional
            Exact gene names to draw; all genes are drawn when empty.
        """
        self.data = data
        self.color = color
        self.label = label
        self.feature_height = feature_height
        self.chevron_height = chevron_height
        self.chevron_width = chevron_width
        self.chevron_spacing = chevron_spacing
        self.truncation_size = truncation_size
        self.utr_endcap_width = utr_endcap_width
        self.gene_name_fontsize = gene_name_fontsize
        self.genes_include = set(genes_include or [])

    @staticmethod
    def _arrange_genes(gene_data_list: Iterable[dict]) -> list[list[str]]:
        """Assign genes to levels so that no two genes on a level overlap.

        Parameters
        ----------
        gene_data_list : iterable of dict
            Gene records with ``ID``, ``start`` and ``end``.

        Returns
        -------
        list of list of str
            Gene IDs per level, longest genes placed first. Always at least one level.
        """
        gene_data_list = sorted(
            gene_data_list, key=lambda gene: gene["end"] - gene["start"], reverse=True
        )
        display_levels = [intervaltree.IntervalTree()]
        for gene_data in gene_data_list:
            level_idx = 0
            while display_levels[level_idx].overlaps(gene_data["start"], gene_data["end"]):
                level_idx += 1
                if level_idx == len(display_levels):
                    display_levels.append(intervaltree.IntervalTree())
            display_levels[level_idx].addi(gene_data["start"], gene_data["end"], data=gene_data)

        return [[gene_interval.data["ID"] for gene_interval in level] for level in display_levels]

    def _sizes_in_data_units(
        self, ax: matplotlib.axes.Axes, ws: int, we: int, fig_width: float, row_height: float
    ) -> dict[str, float]:
        """Convert the inch-valued glyph sizes to data units for this window and row."""
        x_per_inch = (we - ws) / fig_width
        y_per_inch = (ax.get_ylim()[1] - ax.get_ylim()[0]) / row_height
        return {
            "chevron_spacing": self.chevron_spacing * x_per_inch,
            "chevron_width": self.chevron_width * x_per_inch,
            "truncation_width": self.truncation_size * x_per_inch,
            "utr_endcap_width": self.utr_endcap_width * x_per_inch,
            "feature_height": self.feature_height * y_per_inch,
            "chevron_height": self.chevron_height * y_per_inch,
            "truncation_height": self.truncation_size * y_per_inch,
        }

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        """Draw every gene overlapping the window, each with its own components only."""
        genes, transcripts, components, _ = self.data.query(chrom, ws, we)
        if self.genes_include:
            genes = {
                gene_id: gene
                for gene_id, gene in genes.items()
                if gene["Name"] in self.genes_include
            }
        logger.info("Drawing %d genes in %s:%d-%d", len(genes), chrom, ws, we)

        gene_display_levels = self._arrange_genes(genes.values())
        ax.set_ylim((-0.5, len(gene_display_levels) - 0.5))
        sizes = self._sizes_in_data_units(ax, ws, we, fig_width, row_height)

        for level, level_gene_ids in enumerate(gene_display_levels):
            for gene_id in level_gene_ids:
                gene = genes[gene_id]
                self._draw_gene_body(ax, gene, level, ws, we, sizes)
                component_ids = {
                    component_id
                    for transcript_id in gene["transcripts"]
                    for component_id in transcripts[transcript_id]["components"]
                }
                for component_id in sorted(component_ids):
                    self._draw_component(ax, components[component_id], gene["strand"], level, sizes)

        ax.set_yticks([])
        ax.set_ylabel(self.label)

    def _draw_gene_body(
        self,
        ax: matplotlib.axes.Axes,
        gene: dict,
        level: int,
        ws: int,
        we: int,
        sizes: dict[str, float],
    ) -> None:
        """Draw the backbone, name, strand chevrons and truncation marks of one gene."""
        left_truncated = gene["start"] < ws
        right_truncated = gene["end"] > we
        visible_start = max(gene["start"], ws) + left_truncated * sizes["truncation_width"] * 2
        visible_end = min(gene["end"], we) - right_truncated * sizes["truncation_width"] * 2

        ax.plot((visible_start, visible_end), (level, level), color=self.color)
        ax.text(
            x=(visible_start + visible_end) / 2,
            y=level + sizes["feature_height"] * 1.5,
            s=gene["Name"],
            ha="center",
            fontsize=self.gene_name_fontsize,
        )

        span = visible_end - visible_start
        num_chevrons = int(max(span / sizes["chevron_spacing"], 1))
        chevron_remainder = span - (num_chevrons - 1) * sizes["chevron_spacing"]
        chevron_dx = -sizes["chevron_width"] if gene["strand"] == "+" else sizes["chevron_width"]
        for chevron_idx in range(num_chevrons):
            chevron_x = (
                visible_start + chevron_idx * sizes["chevron_spacing"] + chevron_remainder / 2
            )
            for chevron_dy in (sizes["chevron_height"], -sizes["chevron_height"]):
                ax.plot(
                    (chevron_x, chevron_x + chevron_dx),
                    (level, level + chevron_dy),
                    color=self.color,
                )

        if left_truncated:
            self._draw_truncation_marks(ax, level, edge_x=ws, direction=1, sizes=sizes)
        if right_truncated:
            self._draw_truncation_marks(ax, level, edge_x=we, direction=-1, sizes=sizes)

    @staticmethod
    def _draw_truncation_marks(
        ax: matplotlib.axes.Axes,
        level: int,
        edge_x: float,
        direction: int,
        sizes: dict[str, float],
    ) -> None:
        """Draw two arrowheads pointing off a window edge (direction +1 = left edge)."""
        width, height = sizes["truncation_width"], sizes["truncation_height"]
        tip_x = edge_x + direction
        base_x = edge_x + direction * (width + 1)
        x_points = numpy.array([tip_x, base_x, base_x])
        y_points = [level, level - height, level + height]
        for shift in (0, direction * width):
            ax.add_patch(
                matplotlib.patches.Polygon(
                    numpy.vstack([x_points + shift, y_points]).T,
                    edgecolor="k",
                    facecolor="w",
                    fill=True,
                    transform=ax.transData,
                    zorder=3,
                )
            )

    def _draw_component(
        self,
        ax: matplotlib.axes.Axes,
        component: dict,
        strand: str,
        level: int,
        sizes: dict[str, float],
    ) -> None:
        """Draw a CDS as a full-height box, or a UTR as a half-height box with an endcap.

        The endcap sits on the UTR side next to the CDS: the right side for a + strand 5' UTR
        or a - strand 3' UTR, the left side otherwise.
        """
        start, end = component["start"], component["end"]
        feature_height = sizes["feature_height"]
        if component["type"] == "CDS":
            boxes = [(start, end - start, feature_height * 2)]
        else:
            cap_width = min(sizes["utr_endcap_width"], end - start)
            cap_on_right = (component["type"] == "five_prime_UTR") == (strand == "+")
            body_start = start if cap_on_right else start + cap_width
            cap_start = end - cap_width if cap_on_right else start
            boxes = [
                (body_start, end - start - cap_width, feature_height),
                (cap_start, cap_width, feature_height * 2),
            ]
        for box_start, box_width, box_height in boxes:
            ax.add_patch(
                matplotlib.patches.Rectangle(
                    xy=(box_start, level - box_height / 2),
                    width=box_width,
                    height=box_height,
                    facecolor=self.color,
                )
            )


class HicPlot(_BrowserSubPlot):
    def __init__(
        self,
        data,
        label: str = "",
        vertical_scale: float = 1,
        cmap: str = "YlOrRd",
        label_rotation: float = 0,
        transform: Callable[[pandas.DataFrame], pandas.DataFrame] = lambda x: x**2,
        max_masked_diag: int = 2,
    ):
        """Draw a Hi-C contact matrix rotated 45 degrees, upper triangle only.

        Parameters
        ----------
        data : _MatrixData
            Matrix source with a ``bin_size`` attribute, e.g. ``HicDataDict``.
        label : str
            y-axis label.
        vertical_scale : float
            Reserved; not used yet.
        cmap : str
            Colormap.
        label_rotation : float
            y-label rotation in degrees.
        transform : callable
            Applied to the rotated matrix before drawing.
        max_masked_diag : int
            Cells with ``|row - column| < max_masked_diag`` are set to 0.
        """
        self.data = data
        self.cmap = cmap
        self.transform = transform
        self.max_masked_diag = max_masked_diag
        self.label = label
        self.label_rotation = label_rotation
        self.vertical_scale = vertical_scale

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        """Draw the matrix for the window and return the drawn (transformed) data.

        Returns
        -------
        pandas.DataFrame
            The rotated, trimmed and transformed matrix passed to ``imshow``.
        """
        visible_start_bin = utilities.roundto(ws, self.data.bin_size)
        visible_end_bin = utilities.roundto(we, self.data.bin_size)
        visible_span = visible_end_bin - visible_start_bin

        queried = self.data.query(
            chrom, visible_start_bin - visible_span // 2, visible_end_bin + visible_span // 2
        )
        matrix = queried.to_numpy(dtype=float, copy=True)
        row_idx, col_idx = numpy.indices(matrix.shape)
        matrix[numpy.abs(row_idx - col_idx) < self.max_masked_diag] = 0

        plot_data = pandas.DataFrame(
            ndimage.rotate(matrix, 45, reshape=False),
            index=queried.index,
            columns=queried.columns,
        )
        plot_data = plot_data.loc[
            visible_start_bin:visible_end_bin, visible_start_bin:visible_end_bin
        ]
        # Rotated 45 degrees, the upper half of the rows holds the upper triangle.
        plot_data = self.transform(plot_data.iloc[: -plot_data.shape[0] // 2, :])

        ax.set_ylim(0, plot_data.shape[1])
        ax.imshow(plot_data, cmap=self.cmap, aspect="auto", extent=(ws, we, 0, plot_data.shape[1]))
        ax.set_xticks([])
        ax.set_ylabel(self.label, rotation=self.label_rotation, labelpad=DEFAULT_YLABEL_PAD)
        return plot_data


def compute_figure_height(row_heights: Sequence[float], ax_spacing: float) -> float:
    """Return the figure height in inches for the given rows and spacing.

    Parameters
    ----------
    row_heights : sequence of float
        Height of each row, in inches.
    ax_spacing : float
        Gap below each row, as a fraction of the summed row heights.

    Returns
    -------
    float
        Figure height in inches.
    """
    return float(numpy.sum(row_heights)) * (1 + ax_spacing * len(row_heights))


def compute_ax_row_positions(
    row_heights: Sequence[float], ax_spacing: float = 0.1
) -> tuple[list[float], list[float]]:
    """Lay out rows top to bottom in figure coordinates.

    Parameters
    ----------
    row_heights : sequence of float
        Height of each row, in inches.
    ax_spacing : float
        Gap below each row, as a fraction of the summed row heights. Converted to figure
        coordinates with the same figure height that ``compute_figure_height`` returns.

    Returns
    -------
    tuple of (list of float, list of float)
        Bottom coordinate and height of each row, in figure coordinates in [0, 1], suitable
        for ``Figure.add_axes``.
    """
    fig_height = compute_figure_height(row_heights=row_heights, ax_spacing=ax_spacing)
    gap_height = ax_spacing * float(numpy.sum(row_heights)) / fig_height
    heights = [row_height / fig_height for row_height in row_heights]
    bottoms = []
    current_top = 1.0
    for height in heights:
        bottoms.append(current_top - height)
        current_top -= height + gap_height
    return bottoms, heights


class GenomeBrowser:
    def __init__(self, subplot_objects: Sequence[Sequence[_BrowserSubPlot]]):
        """Stack genomic tracks into one figure.

        Parameters
        ----------
        subplot_objects : sequence of sequence of _BrowserSubPlot
            One inner sequence per row, top to bottom. Every track in a row draws on the
            same Axes.
        """
        self.subplot_objects = subplot_objects

    def visualize(
        self,
        chrom: str,
        start: int,
        end: int,
        fig_width: float = 12,
        row_heights: float | Sequence[float] = 1,
        ax_spacing: float = 0.05,
        num_xticks: int = 10,
        style: dict | None = None,
    ) -> matplotlib.figure.Figure:
        """Draw every row for one genomic window and return the figure.

        Parameters
        ----------
        chrom : str
            Chromosome, UCSC style (``chr1``).
        start, end : int
            Window, 0-based half-open.
        fig_width : float
            Figure width in inches.
        row_heights : float or sequence of float
            Row height in inches, one value for all rows or one per row.
        ax_spacing : float
            Gap below each row, as a fraction of the summed row heights.
        num_xticks : int
            Approximate number of x ticks on the bottom row.
        style : dict, optional
            matplotlib rcParams that override ``DEFAULT_STYLE`` for this figure only.

        Returns
        -------
        matplotlib.figure.Figure
            A new figure with one Axes per row.

        Raises
        ------
        ValueError
            If ``end <= start`` or the number of row heights does not match the rows.
        """
        start, end = int(start), int(end)
        if end <= start:
            raise ValueError(f"Window end must be greater than window start; got {start}, {end}")
        if numpy.ndim(row_heights) == 0:
            row_heights = [row_heights] * len(self.subplot_objects)
        if len(row_heights) != len(self.subplot_objects):
            raise ValueError(
                f"Got {len(row_heights)} row heights for {len(self.subplot_objects)} rows"
            )
        logger.info("Drawing %d rows for %s:%d-%d", len(self.subplot_objects), chrom, start, end)

        with matplotlib.rc_context(rc={**DEFAULT_STYLE, **(style or {})}):
            fig = plt.figure(
                figsize=(
                    fig_width,
                    compute_figure_height(row_heights=row_heights, ax_spacing=ax_spacing),
                )
            )
            bottoms, heights = compute_ax_row_positions(
                row_heights=row_heights, ax_spacing=ax_spacing
            )
            for row_idx, row_subplots in enumerate(self.subplot_objects):
                row_ax = fig.add_axes((0, bottoms[row_idx], 1, heights[row_idx]))
                self._set_x_axis(row_ax, chrom, start, end, num_xticks, row_idx)
                # Default limits; tracks may change them.
                row_ax.set_ylim(0, 1)
                row_ax.set_xlim(start, end)
                for subplot in row_subplots:
                    subplot.plot(
                        row_ax,
                        chrom=chrom,
                        ws=start,
                        we=end,
                        fig_width=fig_width,
                        row_height=row_heights[row_idx],
                    )
        logger.info("Done drawing %d rows", len(self.subplot_objects))
        return fig

    def _set_x_axis(
        self,
        ax: matplotlib.axes.Axes,
        chrom: str,
        start: int,
        end: int,
        num_xticks: int,
        row_idx: int,
    ) -> None:
        """Label genomic positions on the bottom row and clear x ticks on the others."""
        if row_idx < len(self.subplot_objects) - 1:
            ax.set_xlabel("")
            ax.set_xticks([])
            return
        span = end - start
        xtick_increment = span / num_xticks
        rounding_increment = 5 * 10 ** numpy.round(numpy.log10(xtick_increment) - 1)
        xtick_increment = utilities.roundto(xtick_increment, rounding_increment)
        num_ticks = int(span / xtick_increment) + 1
        round_start = utilities.roundto(start, rounding_increment)
        ax.set_xticks(numpy.arange(num_ticks) * xtick_increment + round_start)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
        ax.set_xlabel(f"{chrom} position")


def match_ylims(fig: matplotlib.figure.Figure, ax_nums: Sequence[int]) -> None:
    """Give the listed Axes of a figure the same y range, from 0 to the largest upper limit.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure returned by ``GenomeBrowser.visualize``.
    ax_nums : sequence of int
        Indices into ``fig.get_axes()``.
    """
    max_extent = max(fig.get_axes()[ax_num].get_ylim()[1] for ax_num in ax_nums)
    for ax_num in ax_nums:
        fig.get_axes()[ax_num].set_ylim((0, max_extent))
