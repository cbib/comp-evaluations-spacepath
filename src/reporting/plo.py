import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.DataFrame({
    "eval_type": ["A"] * 20 + ["B"] * 20,
    "set_type": ["pathway"] * 10 + ["ghost"] * 10 + ["pathway"] * 10 + ["ghost"] * 10,
    "logfoldchanges": range(40)
})

palette_set_type_dark = {"ghost" : "#fc8d62",   "pathway" : "#66c2a5"}
palette_set_type_pastel = {"ghost" : "#FDCDAC",   # pastel orange/coral
                           "pathway": "#B3E2CD"}  # pastel teal



sns.boxplot(
    data=df,
    x="eval_type",
    y="logfoldchanges",
    hue="set_type",
    palette=palette_set_type
)

plt.show()