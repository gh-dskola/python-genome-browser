"""Polygon approximations of elliptical arcs, drawn with ``Axes.plot``.

matplotlib's ``Arc`` patch is not used because it rendered incorrectly when this code was
written. Angles follow the original convention: theta runs from pi/2 (left end of the arc)
through pi (apex) to 3*pi/2 (right end).
"""

import matplotlib.axes
import numpy

DEFAULT_ARC_POINTS = 200
LEFT_END_THETA = numpy.pi / 2
RIGHT_END_THETA = numpy.pi * 1.5
MIDPOINT_THETA = numpy.pi
# Keeps arccos-derived angles strictly inside the arc's end points.
THETA_EPSILON = 0.00001


def compute_half_arc_points(
    center: tuple[float, float],
    a: float,
    b: float,
    theta1: float,
    theta2: float,
    num_points: int = DEFAULT_ARC_POINTS,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Compute points on one quadrant of an ellipse.

    Parameters
    ----------
    center : tuple of float
        (x, y) center of the ellipse.
    a, b : float
        Horizontal and vertical semi-axes.
    theta1, theta2 : float
        First and last angle, within one quadrant.
    num_points : int
        Number of points to return.

    Returns
    -------
    tuple of numpy.ndarray
        x and y coordinates.
    """
    theta = theta1 + (theta2 - theta1) * numpy.arange(num_points) / max(num_points - 1, 1)
    phi = numpy.pi / 2 - numpy.arctan(numpy.tan(theta))
    return center[0] + a * numpy.cos(phi), center[1] + b * numpy.sin(phi)


def draw_arc(
    ax: matplotlib.axes.Axes,
    center: tuple[float, float],
    height: float,
    width: float,
    theta1: float = LEFT_END_THETA,
    theta2: float = RIGHT_END_THETA,
    color: object = "k",
    direction: str = "down",
    num_points: int = DEFAULT_ARC_POINTS,
) -> None:
    """Draw the part of a half-ellipse between two angles.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    center : tuple of float
        (x, baseline y) of the arc.
    height, width : float
        Arc height and full width in data units.
    theta1, theta2 : float
        Angle range, with ``LEFT_END_THETA <= theta1 <= theta2 <= RIGHT_END_THETA``.
    color : matplotlib color
        Line color.
    direction : str
        ``down`` draws below the baseline, anything else above it.
    num_points : int
        Number of polygon points.
    """
    if not LEFT_END_THETA <= theta1 <= theta2 <= RIGHT_END_THETA:
        raise ValueError(f"Arc angles out of range: theta1={theta1}, theta2={theta2}")

    left_angle_span = min(max(MIDPOINT_THETA - theta1, 0), theta2 - theta1)
    right_angle_span = min(max(theta2 - MIDPOINT_THETA, 0), theta2 - theta1)
    left_points = int(num_points * left_angle_span / (left_angle_span + right_angle_span))
    right_points = num_points - left_points

    segments = []
    if left_points:
        segments.append(
            compute_half_arc_points(
                (center[0], 0), width / 2, height, theta1, theta1 + left_angle_span, left_points
            )
        )
    if right_points:
        segments.append(
            compute_half_arc_points(
                (center[0], 0), width / 2, height, theta2 - right_angle_span, theta2, right_points
            )
        )
    x_coords = numpy.concatenate([segment[0] for segment in segments])
    y_coords = numpy.concatenate([segment[1] for segment in segments])
    if direction == "down":
        y_coords = -y_coords

    ax.plot(x_coords, y_coords + center[1], color=color)


def draw_visible_arc(
    ax: matplotlib.axes.Axes,
    center: tuple[float, float],
    height: float,
    width: float,
    ws: float,
    we: float,
    color: object = "k",
    direction: str = "down",
    num_points: int = DEFAULT_ARC_POINTS,
) -> None:
    """Draw a half-ellipse clipped to the x range [ws, we].

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    center : tuple of float
        (x, baseline y) of the arc.
    height, width : float
        Arc height and full width in data units.
    ws, we : float
        Visible x range.
    color : matplotlib color
        Line color.
    direction : str
        ``down`` draws below the baseline, anything else above it.
    num_points : int
        Number of polygon points.
    """

    def theta_at(x: float) -> float:
        """Return the arc angle whose point has x coordinate ``x``."""
        phi = numpy.arccos((x - center[0]) / (width / 2))
        return numpy.arctan(1 / numpy.tan(phi)) + numpy.pi

    theta_start = theta_at(ws) if ws > center[0] - width / 2 else LEFT_END_THETA + THETA_EPSILON
    theta_end = theta_at(we) if we < center[0] + width / 2 else RIGHT_END_THETA - THETA_EPSILON

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
    ax: matplotlib.axes.Axes,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    left_bin_center: float,
    right_bin_center: float,
    color: object = "k",
    baseline: float = 0.0,
    vertical_scaling_factor: float = 1,
    direction: str = "down",
    num_points: int = DEFAULT_ARC_POINTS,
) -> None:
    """Draw one interaction as an arc whose height grows with its span.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    xlim, ylim : tuple of float
        Visible x and y ranges.
    left_bin_center, right_bin_center : float
        x positions of the two anchors.
    color : matplotlib color
        Line color.
    baseline : float
        y position of the anchors.
    vertical_scaling_factor : float
        Multiplier on the arc height.
    direction : str
        ``down`` draws below the baseline, anything else above it.
    num_points : int
        Number of polygon points.
    """
    arc_width = right_bin_center - left_bin_center
    vertical_span = baseline - ylim[0] if direction == "down" else ylim[1] - baseline
    arc_height = vertical_span * (arc_width / (xlim[1] - xlim[0])) * vertical_scaling_factor

    draw_visible_arc(
        ax=ax,
        center=((left_bin_center + right_bin_center) / 2, baseline),
        height=arc_height,
        width=arc_width,
        ws=xlim[0],
        we=xlim[1],
        color=color,
        direction=direction,
        num_points=num_points,
    )
