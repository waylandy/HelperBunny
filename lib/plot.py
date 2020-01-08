import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_ss(ss, x, ypos=20, height=2, ax=None):
    ypos -= height/2
    for s, p in zip(ss, x):
        if s == '-':
            continue
        elif s in 'EB':  # blue
            pat=patches.Rectangle(width=1,height=height,fc='#0000ff',lw=0,xy=(p-0.5, ypos))
        elif s in 'GHI': # red
            pat=patches.Rectangle(width=1,height=height,fc='#ff0000',lw=0,xy=(p-0.5, ypos))
        elif s in 'STC': # black
            pat=patches.Rectangle(width=1,height=height/4,fc='#000000',lw=0, xy=(p-0.5,ypos+(height*(3/8))))
        if ax == None:
            plt.gca().add_patch(pat)
        else:
            ax.add_patch(pat)
