import intervaltree
import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas
from scipy import ndimage

from . import utilities

DEFAULT_ARC_POINTS = 200
DEFAULT_YLABEL_PAD = 50
CHROMOSOME_DIALECT = "ucsc"


# ToDo: Move this stuff to a separate module so we can delete wholesale when the Ellipse class gets fixed.
def compute_half_arc_points(center, a, b, theta1, theta2, num_points=DEFAULT_ARC_POINTS):
    """
    Computes the coordinates for component points of a polygonal approximation to
    an ellipse for a single quadrant.
    """
    # ToDo: Add input validation to make sure we stay within a single quadrant.
    x_coords = numpy.empty(num_points)
    y_coords = numpy.empty(num_points)

    for i in range(0, num_points):
        theta = (theta2 - theta1) * (i / max(num_points - 1, 1)) + theta1
        fi = numpy.pi / 2 - numpy.arctan(numpy.tan(theta))
        x = center[0] + a * numpy.cos(fi)
        y = center[1] + b * numpy.sin(fi)
        x_coords[i] = x
        y_coords[i] = y

    return x_coords, y_coords


def draw_arc(
    ax,
    center,
    height,
    width,
    theta1=0,
    theta2=numpy.pi,
    color="k",
    direction="down",
    num_points=DEFAULT_ARC_POINTS,
):
    """
    Since Matplotlib's Arc Patches are broken at the moment, we draw arcs using the ax.plot() method
    instead.
    """
    LEFT_END_THETA = numpy.pi / 2
    RIGHT_END_THETA = numpy.pi * 1.5
    MIDPOINT_THETA = numpy.pi

    vertical_baseline = center[1]

    assert LEFT_END_THETA <= theta1 <= theta2 <= RIGHT_END_THETA

    b = height
    a = width / 2

    # determine how to allocate points
    left_angle_span = min(max(MIDPOINT_THETA - theta1, 0), theta2 - theta1)
    right_angle_span = min(max(theta2 - MIDPOINT_THETA, 0), theta2 - theta1)
    total_angle_span = left_angle_span + right_angle_span
    left_points = int(num_points * left_angle_span / total_angle_span)
    right_points = num_points - left_points

    x_coords = numpy.empty(num_points)
    y_coords = numpy.empty(num_points)

    if left_points:
        # plot upper left quadrant
        left_theta2 = theta1 + left_angle_span
        x, y = compute_half_arc_points(
            center=(center[0], 0),
            a=a,
            b=b,
            theta1=theta1,
            theta2=left_theta2,
            num_points=left_points,
        )
        x_coords[:left_points] = x[:]
        y_coords[:left_points] = y[:]
    if right_points:
        # plot upper right quadrant
        right_theta1 = theta2 - right_angle_span
        x, y = compute_half_arc_points(
            center=(center[0], 0),
            a=a,
            b=b,
            theta1=right_theta1,
            theta2=theta2,
            num_points=right_points,
        )
        x_coords[left_points:] = x[:]
        y_coords[left_points:] = y[:]

    if direction == "down":
        y_coords = -y_coords

    y_coords += vertical_baseline

    ax.plot(x_coords, y_coords, color=color)


def draw_visible_arc(
    ax, center, height, width, ws, we, color="k", direction="down", num_points=DEFAULT_ARC_POINTS
):
    """
    Draws a 180 degree elliptical arc truncated by an interval of x coordinates.
    Does not truncate based on y coordinates
    """
    # ToDo: Subtract 1 pi from all coordinates
    LEFT_END_THETA = numpy.pi / 2 + 0.00001
    RIGHT_END_THETA = numpy.pi * 1.5 - 0.00001

    def infer_theta_cutoff(x, arc_center, arc_width):
        a = arc_width / 2
        fi = numpy.arccos((x - arc_center[0]) / a)
        theta = numpy.arctan(1 / numpy.tan(fi)) + numpy.pi
        return theta

    if ws > center[0] - width / 2:
        theta_start = infer_theta_cutoff(x=ws, arc_center=center, arc_width=width)
    else:
        theta_start = LEFT_END_THETA

    if we < center[0] + width / 2:
        theta_end = infer_theta_cutoff(x=we, arc_center=center, arc_width=width)
    else:
        theta_end = RIGHT_END_THETA

    draw_arc(
        ax=ax,
        center=center,
        height=height,
        width=width,
        theta1=theta_start,
        theta2=theta_end,
        color=color,
        direction=direction,
        num_points=num_points,
    )


def draw_arc_interaction(
    ax,
    xlim,
    ylim,
    left_bin_center,
    right_bin_center,
    color="k",
    baseline=0.0,
    vertical_scaling_factor=1,
    direction="down",
    num_points=DEFAULT_ARC_POINTS,
):
    """ """
    arc_width = right_bin_center - left_bin_center

    if direction == "down":
        vertical_span = baseline - ylim[0]
    else:
        vertical_span = ylim[1] - baseline

    arc_height = vertical_span * (arc_width / (xlim[1] - xlim[0])) * vertical_scaling_factor

    draw_visible_arc(
        ax=ax,
        center=((left_bin_center + right_bin_center) / 2, baseline),
        height=arc_height,
        ws=xlim[0],
        we=xlim[1],
        width=arc_width,
        color=color,
        direction=direction,
        num_points=num_points,
    )


class _BrowserSubPlot:
    def __init__(self):
        # self.chrom = None
        # self.ws = None
        # self.we = None
        # self.fig_width = None
        # self.row_height = None
        pass

    # def set_globals(self, chrom, ws, we, fig_width=64, row_height=4):
    #     self.chrom = chrom
    #     self.ws = ws
    #     self.we = we
    #     self.fig_width = fig_width
    #     self.row_height = row_height

    # @property
    # def aspect_ratio(self):
    #     return self.fig_width / self.row_height

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        print("Stub method -- must be overridden by inheritors")


class InteractionPlot(_BrowserSubPlot):
    def __init__(
        self,
        interaction_df,
        bin_size,
        arc_color=(0.7, 0.3, 0.6),
        direction="down",
        baseline=None,
        vertical_scaling_factor=1,
        thickness_column=None,
        show_bin_centers=True,
        label="Plac-seq",
    ):
        super(InteractionPlot, self).__init__()

        self.interaction_df = interaction_df
        self.bin_size = bin_size
        if baseline is None:
            if direction == "down":
                baseline = 1
            else:
                baseline = 0
        self.baseline = baseline
        self.vertical_scaling_factor = vertical_scaling_factor
        self.arc_color = arc_color
        self.label = label
        self.direction = direction
        self.show_bin_centers = show_bin_centers
        self.thickness_column = thickness_column

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        # Filter the interaction DataFrame to interactions with at least one anchor point within the visible window.
        visible_interactions = self.interaction_df.loc[self.interaction_df["chr1"] == chrom]
        left_bin_midpoints = (visible_interactions["end1"] + visible_interactions["start1"]) / 2
        right_bin_midpoints = (visible_interactions["end2"] + visible_interactions["start2"]) / 2
        left_visible = (left_bin_midpoints >= ws) & (left_bin_midpoints <= we)
        right_visible = (right_bin_midpoints >= ws) & (right_bin_midpoints <= we)
        visible_interactions = visible_interactions.loc[left_visible | right_visible]

        original_ylim = ax.get_ylim()
        ax.set_xlim(ws, we)  # ToDo: Standardize this behavior across all subplot classes

        for interaction_id in visible_interactions.index:
            draw_arc_interaction(
                ax,
                left_bin_center=left_bin_midpoints.loc[interaction_id],
                right_bin_center=right_bin_midpoints.loc[interaction_id],
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
            leftmost_tick = (
                numpy.ceil((ws - self.bin_size / 2) / self.bin_size) * self.bin_size
                + self.bin_size / 2
            )
            rightmost_tick = (
                numpy.floor((we - self.bin_size / 2) / self.bin_size + 1) * self.bin_size
                + self.bin_size / 2
            )

            ax.set_xticks(numpy.arange(leftmost_tick, rightmost_tick, self.bin_size))
            ax.set_xticklabels([])

            if self.direction == "down":
                ax.xaxis.set_ticks_position("top")

        ax.set_yticks([])


class BedPlot(_BrowserSubPlot):
    DEFAULT_PATCH_KWARGS = {"linewidth": 1, "edgecolor": "k"}

    def __init__(
        self,
        interval_data,
        label="",
        cmap="RdBu_r",
        baseline=0.5,
        color="k",
        color_by="",
        display_value="",
        patch_height=0.5,
        pad_fraction=0.1,
        patch_kwargs=None,
    ):
        super().__init__()

        self.interval_data = interval_data
        self.color = color
        self.pad_fraction = pad_fraction

        columns = interval_data.data.columns
        self.color_by = color_by
        if self.color_by:
            if self.color_by not in columns:
                raise ValueError(f"Color-by column {self.color_by!r} is not in the interval data")
            extent = numpy.abs(interval_data.data[self.color_by]).max()
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
        utilities.add_label(ax=ax, tick=self.baseline, tick_label=self.label, axis="y")

        visible_intervals = self.interval_data.query(chrom, ws, we)
        for _, this_interval in visible_intervals.iterrows():
            start_loc = this_interval["chromStart"]
            end_loc = this_interval["chromEnd"]
            if end_loc <= start_loc:
                raise ValueError(f"Interval end {end_loc} is not greater than start {start_loc}")

            if self.color_by:
                interval_color = self.color_mapper.to_rgba(this_interval[self.color_by])
            else:
                interval_color = self.color

            rec = matplotlib.patches.Rectangle(
                xy=(start_loc, self.baseline - self.patch_height / 2),
                width=end_loc - start_loc,
                height=self.patch_height,
                facecolor=interval_color,
                **self.patch_kwargs,
            )
            ax.add_patch(rec)
            if self.display_value:
                ax.text(
                    x=(start_loc + end_loc) / 2,
                    y=self.baseline,
                    s=f"{this_interval[self.display_value]:>0.2}",
                    ha="center",
                )

        for edge in (self.baseline + self.patch_height / 2, self.baseline - self.patch_height / 2):
            utilities.adjust_limits(
                ax=ax, new_position=edge, axis="y", padding_fraction=self.pad_fraction
            )


class WigPlot(_BrowserSubPlot):
    # ToDo: Add support for stranded data
    def __init__(
        self,
        genomic_vector_data,
        label=None,
        color=None,
        solid=True,
        alpha=1.0,
        center_vector=False,
        scale_vector_to_plot=False,
        label_rotation=0,
    ):
        super().__init__()
        self.data = genomic_vector_data
        self.color = color
        self.solid = solid
        self.alpha = alpha
        self.center = center_vector
        self.scale_vector_to_plot = scale_vector_to_plot
        self.label = label
        self.label_rotation = label_rotation

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        ylim = ax.get_ylim()
        vert_span = ylim[1] - ylim[0]
        vert_center = vert_span / 2 + ylim[0]

        # astype returns a copy, so scaling below never writes into the data source.
        this_plot_vector = self.data[chrom].loc[ws:we].astype(float)

        if self.scale_vector_to_plot:
            this_plot_vector = (
                this_plot_vector / (this_plot_vector.max() - this_plot_vector.min()) * vert_span
            )
        if self.center:
            this_plot_vector = this_plot_vector - this_plot_vector.mean() + vert_center

        this_plot_vector = this_plot_vector.loc[
            (this_plot_vector.index >= ws) & (this_plot_vector.index < we)
        ]
        this_plot_vector.name = self.label

        if self.solid:
            ax.fill_between(
                x=this_plot_vector.index,
                y1=this_plot_vector,
                color=self.color,
                alpha=self.alpha,
                label=self.label,
            )
        else:
            ax.plot(
                this_plot_vector.index,
                this_plot_vector,
                color=self.color,
                alpha=self.alpha,
                label=self.label,
            )

        ax.autoscale(enable=True, axis="y")

        if self.label:
            ax.set_ylabel(self.label, rotation=self.label_rotation, labelpad=DEFAULT_YLABEL_PAD)


class GeneModelPlot(_BrowserSubPlot):
    def __init__(
        self,
        gene_annotation_data,
        label="Genes",
        color="k",
        feature_height=0.12,
        chevron_height=0.05,
        chevron_width=0.04,
        chevron_spacing=0.10,
        truncation_size=0.10,
        utr_endcap_width=0.04,
        gene_name_fontsize=8,
        genes_include=None,
    ):
        super().__init__()

        self.gene_annotation_data = gene_annotation_data
        self.color = color
        self.label = label
        self.feature_height = feature_height  # in inches
        self.chevron_height = chevron_height  # in inches
        self.chevron_width = chevron_width  # in inches
        self.chevron_spacing = chevron_spacing  # in inches
        self.truncation_size = truncation_size  # in inches
        self.utr_endcap_width = utr_endcap_width  # in inches
        self.gene_name_fontsize = gene_name_fontsize
        self.genes_include = set(genes_include or [])

    @staticmethod
    def _arrange_genes(gene_data_list):
        """
        Given an iterable of gene data dictionaries,
        returns a list of lists of gene names that
        should be displayed at various levels.
        """
        gene_data_list = sorted(gene_data_list, key=lambda x: x["end"] - x["start"], reverse=True)

        display_levels = [intervaltree.IntervalTree()]
        for gene_data in gene_data_list:
            level_idx = 0
            while display_levels[level_idx].overlaps(gene_data["start"], gene_data["end"]):
                level_idx += 1
                if level_idx == len(display_levels):
                    display_levels.append(intervaltree.IntervalTree())
            display_levels[level_idx].addi(gene_data["start"], gene_data["end"], data=gene_data)

        return [[gene_interval.data["ID"] for gene_interval in level] for level in display_levels]

    def _sizes_in_data_units(self, ax, ws, we, fig_width, row_height):
        """Convert the inch-valued glyph sizes to data coordinates for this window and row."""
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
        genes, transcripts, components, _ = self.gene_annotation_data.query(chrom, ws, we)
        if self.genes_include:
            genes = {
                gene_id: gene for gene_id, gene in genes.items()
                if gene["Name"] in self.genes_include
            }

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

    def _draw_gene_body(self, ax, gene, level, ws, we, sizes):
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
    def _draw_truncation_marks(ax, level, edge_x, direction, sizes):
        """Draw two arrowheads pointing off the window edge (direction +1 = left edge)."""
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

    def _draw_component(self, ax, component, strand, level, sizes):
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


def compute_figure_height(row_heights: list[float], ax_spacing: float) -> float:
    """Return the figure height in inches for the given rows and spacing.

    Parameters
    ----------
    row_heights : list of float
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
    row_heights: list[float], ax_spacing: float = 0.1
) -> tuple[list[float], list[float]]:
    """Lay out rows top to bottom in figure coordinates.

    Parameters
    ----------
    row_heights : list of float
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
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans", "Bitstream Vera Sans", "sans-serif"],
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


class GenomeBrowser:
    def __init__(self, subplot_objects):
        """Stack genomic tracks into one figure.

        Parameters
        ----------
        subplot_objects : list of list of _BrowserSubPlot
            One inner list per row, top to bottom. Every track in a row draws on the same Axes.
        """
        self.subplot_objects = subplot_objects

    def visualize(
        self,
        chrom,
        start,
        end,
        fig_width=12,
        row_heights=1,
        ax_spacing=0.05,
        num_xticks=10,
        style=None,
    ):
        """Draw every row for one genomic window and return the figure.

        Parameters
        ----------
        chrom : str
            Chromosome, UCSC style (``chr1``).
        start, end : int
            Window, 0-based half-open.
        fig_width : float
            Figure width in inches.
        row_heights : float or list of float
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
                row_ax = fig.add_axes([0, bottoms[row_idx], 1, heights[row_idx]])
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
        return fig

    def _set_x_axis(self, ax, chrom, start, end, num_xticks, row_idx):
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


class HicPlot(_BrowserSubPlot):
    def __init__(
        self,
        data,
        label="",
        vertical_scale=1,
        cmap="YlOrRd",
        label_rotation=0,
        transform=lambda x: x**2,
        max_masked_diag=2,
    ):
        super().__init__()
        self.data = data
        self.cmap = cmap
        self.transform = transform
        self.max_masked_diag = max_masked_diag
        self.label = label
        self.label_rotation = label_rotation
        self.vertical_scale = vertical_scale

    def plot(self, ax, chrom, ws, we, fig_width, row_height):
        visible_start_bin = utilities.roundto(ws, self.data.bin_size)
        visible_end_bin = utilities.roundto(we, self.data.bin_size)
        visible_span = visible_end_bin - visible_start_bin

        data_start_bin = visible_start_bin - visible_span // 2
        data_end_bin = visible_end_bin + visible_span // 2

        queried = self.data.query(chrom, data_start_bin, data_end_bin)
        matrix = queried.to_numpy(dtype=float, copy=True)
        row_idx, col_idx = numpy.indices(matrix.shape)
        matrix[numpy.abs(row_idx - col_idx) < self.max_masked_diag] = 0

        plot_data = pandas.DataFrame(
            ndimage.rotate(matrix, 45, reshape=False),
            index=queried.index,
            columns=queried.columns,
        )
        plot_data = plot_data.loc[visible_start_bin:visible_end_bin, visible_start_bin:visible_end_bin]
        # Rotated 45 degrees, the upper half of the rows holds the upper triangle.
        plot_data = plot_data.iloc[: -plot_data.shape[0] // 2, :]
        plot_data = self.transform(plot_data)

        ax.set_ylim(0, plot_data.shape[1])
        ax.imshow(
            plot_data, cmap=self.cmap, aspect="auto", extent=(ws, we, 0, plot_data.shape[1])
        )
        ax.set_xticks([])
        ax.set_ylabel(self.label, rotation=self.label_rotation, labelpad=DEFAULT_YLABEL_PAD)
        return plot_data


def match_ylims(fig, ax_nums):
    """
    Will make the upper ylim of each of the numbered axes of :param fig: listed in
    :param ax_nums: equal to the maximum found in any of the numbered axes.
    """
    max_extent = max([fig.get_axes()[ax_num].get_ylim()[1] for ax_num in ax_nums])
    for ax_num in ax_nums:
        fig.get_axes()[ax_num].set_ylim((0, max_extent))
