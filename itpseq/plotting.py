import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.collections import LineCollection

__all__ = ['itoeprint_plot']


def itoeprint_plot(
    dataset=None,
    plot='bands',
    norm='sum',
    norm_range=(21, 51),
    exposure=1,
    limit=(0, 100),
    show_range=True,
    ax=None,
):
    if not ax:
        ax = plt.subplot()

    df = dataset.toeprint_df

    start, end = norm_range
    low_limit, high_limit = limit

    df = df.div(df.loc[end:start].agg(norm) if norm else df.max().max() / 20)

    xticks_kwargs = dict(
        rotation='vertical' if df.shape[1] > 10 else 'horizontal', ha='center'
    )

    if plot == 'bands':
        for i, c in enumerate(df):
            s = df[c].dropna()
            points = [[(i - 0.4, y), (i + 0.4, y)] for y in s.index]
            # ax.add_collection(LineCollection(points, color='k', lw=s*10, alpha=0.3))
            ax.add_collection(
                LineCollection(
                    points,
                    color='k',
                    lw=s * 0.2 * exposure,
                    path_effects=[path_effects.Stroke(capstyle='round')],
                    # path_effects=[path_effects.Stroke(capstyle='butt')],
                )
            )

        ax.set_xlim(-0.5, df.shape[1] - 0.5)
        ax.set_ylim(low_limit, high_limit)
        ax.set_xticks(np.arange(df.shape[1]), df.columns, **xticks_kwargs)

        ax.spines[['left', 'right', 'top', 'bottom']].set_visible(False)

        if show_range:
            ax.axhline(start, color='red', alpha=0.5)
            ax.axhline(end, color='red', alpha=0.5)
            ax.axhspan(low_limit, start, color='red', alpha=0.1, zorder=2)
            ax.axhspan(end, high_limit, color='red', alpha=0.1, zorder=2)

    elif plot == 'shades':
        sns.heatmap(
            df.pipe(
                lambda x: x.reindex(
                    range(high_limit or df.index.max(), low_limit, -1)
                )
            ),
            cmap='Greys',
            ax=ax,
            cbar=False,
        )
        for i in range(0, df.shape[1] + 1):
            ax.axvline(i, c='w', lw=70 / (df.shape[1] + 1), zorder=1)

        ax.set_xlim(0, i)
        ax.set_xticks(
            np.arange(df.shape[1]) + 0.5, df.columns, **xticks_kwargs
        )

        if show_range:
            ax.axhline(high_limit - start, color='red', alpha=0.5)
            ax.axhline(high_limit - end, color='red', alpha=0.5)
            ax.axhspan(
                high_limit,
                high_limit - start,
                color='red',
                alpha=0.1,
                zorder=1,
            )
            ax.axhspan(
                high_limit - end, low_limit, color='red', alpha=0.1, zorder=1
            )

    ax.set_ylabel('Inverse-toeprint length')
    ax.tick_params(top=False, labeltop=True, bottom=False)

    return ax
