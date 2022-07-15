import pandas as pd
import numpy as np
from bisect import bisect_right
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import time

# References
# https://www.ihsenergy.ca/support/documentation_ca/WellTest/2019_1/content/html_files/reference_materials/nomenclature.htm
# https://www.ihsenergy.ca/support/documentation_ca/WellTest/2019_1/content/html_files/analysis_types/conventional_test_analyses/derivative_analyses.htm
# https://www.ihsenergy.ca/support/documentation_ca/WellTest/content/html_files/analysis_types/conventional_test_analyses/radial_flow_analysis.htm
# https://www.ihsenergy.ca/support/documentation_ca/WellTest/2019_1/content/html_files/analysis_types/conventional_test_analyses/afterflow_analysis.htm#Summary_of_Equations_for_Afterflow_Derivative_Analysis


# start and end points from pressure data
def get_limits(p_data):
    # define Draggable lines class
    class draggable_lines:
        def __init__(self, ax, kind, XorY, label, color):
            self.ax = ax
            self.c = ax.get_figure().canvas
            self.o = kind
            self.XorY = XorY
            self.label = label

            if kind == "h":
                x = [-1, 1]
                y = [XorY, XorY]

            elif kind == "v":
                x = [XorY, XorY]
                y = [0.00001, 1000000]
            self.line = lines.Line2D(x, y, picker=5, label=label, color=color)
            self.ax.add_line(self.line)
            self.c.draw_idle()
            self.sid = self.c.mpl_connect('pick_event', self.clickonline)

        def clickonline(self, event):
            if event.artist == self.line:
                self.follower = self.c.mpl_connect("motion_notify_event", self.followmouse)
                self.releaser = self.c.mpl_connect("button_press_event", self.releaseonclick)

        def followmouse(self, event):
            if self.o == "h":
                self.line.set_ydata([event.ydata, event.ydata])
            else:
                self.line.set_xdata([event.xdata, event.xdata])
            self.c.draw_idle()

        def releaseonclick(self, event):
            global DD_start, DD_end, BU_end
            if self.o == "h":
                self.XorY = self.line.get_ydata()[0]
            else:
                self.XorY = self.line.get_xdata()[0]
                if self.label == "DD_Start":
                    DD_start = round(self.XorY)
                elif self.label == "DD_end":
                    DD_end = round(self.XorY)
                elif self.label == "BU_end":
                    BU_end = round(self.XorY)

            self.c.mpl_disconnect(self.releaser)
            self.c.mpl_disconnect(self.follower)

    # plot pressure data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(p_data.index, p_data["Press"], linestyle="none", marker="o", markersize=3)
    DD1 = draggable_lines(ax, "v", np.percentile(p_data.index, 30), label="DD_Start", color="green")
    DD2 = draggable_lines(ax, "v", np.percentile(p_data.index, 60), label="DD_end", color="red")
    BU2 = draggable_lines(ax, "v", np.percentile(p_data.index, 80), label="BU_end", color="black")
    plt.suptitle("Pressure Plot")
    plt.title("Move the lines then click to estimate start and end times of each flow period DD/BU")
    plt.ylim([0, max(p_data["Press"]) + 1000])
    plt.legend()
    plt.grid()
    plt.show()
    try:
        # Filter Draw-down and Buildup data
        DD_data = p_data.loc[(p_data.index >= DD_start) & (p_data.index <= DD_end), ["DateTime", "Press"]].copy()
        DD_data.reset_index(inplace=True, drop=True)
        BU_data = p_data.loc[(p_data.index >= DD_end) & (p_data.index <= BU_end), ["DateTime", "Press"]].copy()
        BU_data.reset_index(inplace=True, drop=True)
        return DD_data, BU_data
    except NameError:
        print("Limits were not selected")


# Define prepare data function
def prepare_data(raw_DD_data, raw_BU_data, params, required_test="BU"):
    # add parameters to dictionary
    params["tp"] = (raw_DD_data.loc[len(raw_DD_data) - 1, "DateTime"] - raw_DD_data.loc[0, "DateTime"]).total_seconds() / 3600
    params["BU_duration"] = (raw_BU_data.loc[len(raw_BU_data) - 1, "DateTime"] - raw_BU_data.loc[0, "DateTime"]).total_seconds() / 3600
    params["test_type"] = required_test
    # prepare Buildup data
    if params["test_type"] == "BU":
        raw_data = raw_BU_data.copy()
        raw_data["p"] = raw_data["Press"]
        raw_data["t"] = raw_data["DateTime"].apply(lambda x:  (x - raw_data.loc[0, "DateTime"]).total_seconds() / 3600)
        raw_data["dp"] = raw_data["p"] - raw_data.loc[0, "p"]
        raw_data["te"] = raw_data["t"] * params["tp"] / (raw_data["t"] + params["tp"])
        params["pwf"] = raw_data.loc[0, "p"]
    # prepare Draw-down data
    else:
        raw_data = raw_DD_data.copy()
        raw_data["p"] = raw_data["Press"]
        raw_data["t"] = raw_data["DateTime"].apply(lambda x:  (x - raw_data.loc[0, "DateTime"]).total_seconds() / 3600)
        raw_data["dp"] = params["Pi"] - raw_data["p"]
        raw_data["te"] = raw_data["t"]
        params["pwf"] = np.mean(raw_data["p"])
    return raw_data, params


# define Bourdet derivative function
def calc_der(raw_data, params, req_L=0.1):
    params["L"] = req_L

    # define binary function for pws search
    def BinarySearch(a, x, b=pd.Series([0])):
        i = bisect_right(a, x)
        if len(b) == 1:
            if i:
                return a[i - 1]
            else:
                return np.nan
        else:
            if i:
                return b[i - 1]
            else:
                return np.nan
    # Define start time
    t1 = time.time()
    # prepare Bourdet derivative
    raw_data["X_C"] = np.log(raw_data["te"])
    raw_data["X_L"] = raw_data["X_C"].apply(lambda x: BinarySearch(raw_data["X_C"], x - params["L"]))
    raw_data["X_R"] = raw_data["X_C"].apply(lambda x: BinarySearch(raw_data["X_C"], x + params["L"]))

    raw_data["P_C"] = raw_data["dp"]
    raw_data["P_L"] = raw_data["X_L"].apply(lambda x: BinarySearch(raw_data["X_C"], x, b=raw_data["P_C"]))
    raw_data["P_R"] = raw_data["X_R"].apply(lambda x: BinarySearch(raw_data["X_C"], x, b=raw_data["P_C"]))

    raw_data["t(ddelP/dt)_L"] = (raw_data["P_C"] - raw_data["P_L"]) / (raw_data["X_C"] - raw_data["X_L"])
    raw_data["t(ddelP/dt)_R"] = (raw_data["P_R"] - raw_data["P_C"]) / (raw_data["X_R"] - raw_data["X_C"])
    raw_data["t(ddelP/dt)_R"] = raw_data.apply(lambda x: x["t(ddelP/dt)_L"] if np.isnan(x["t(ddelP/dt)_R"]) else x["t(ddelP/dt)_R"], axis=1)
    raw_data["t(ddelP/dt)_C"] = (((raw_data["X_R"] - raw_data["X_C"]) * raw_data["t(ddelP/dt)_L"]) + ((raw_data["X_C"] - raw_data["X_L"]) * raw_data["t(ddelP/dt)_R"])) / (raw_data["X_R"] - raw_data["X_L"])
    max_te = max(raw_data["t"])
    raw_data["t(ddelP/dt)_C"] = raw_data.apply(lambda x: np.nan if x["t"] > max_te - params["L"] else x["t(ddelP/dt)_C"], axis=1)
    raw_data["derv"] = raw_data["t(ddelP/dt)_C"]
    # Print Execution time
    print('Execution time:', round(time.time() - t1, 2), 'seconds')

    return raw_data.loc[:, ["t", "p", "dp", "derv"]], params


# define Derivative plot analysis function
def derivative_plot_analysis(derv_data, params):
    # define Class for draggable lines
    class draggable_lines:
        def __init__(self, ax, kind, XorY):
            self.ax = ax
            self.c = ax.get_figure().canvas
            self.o = kind
            self.XorY = XorY
            if kind == "h":
                x = [0.000001, 1000000]
                y = [XorY, XorY]
            elif kind == "v":
                x = [XorY, XorY]
                y = [0.000001, 1000000]
            elif kind == "US":
                y = [XorY, pow(10, (np.log10(XorY) + np.log10(1000000) - np.log10(0.0001)))]
                x = [0.0001, 1000000]
            self.line = lines.Line2D(x, y, picker=5, color="black", linestyle="dashed")
            self.ax.add_line(self.line)
            self.c.draw_idle()
            self.sid = self.c.mpl_connect('pick_event', self.clickonline)

        def clickonline(self, event):
            if event.artist == self.line:
                self.follower = self.c.mpl_connect("motion_notify_event", self.followmouse)
                self.releaser = self.c.mpl_connect("button_press_event", self.releaseonclick)

        def followmouse(self, event):
            if self.o == "h":
                self.line.set_ydata([event.ydata, event.ydata])
            elif self.o == "v":
                self.line.set_xdata([event.xdata, event.xdata])
            elif self.o == "US":
                self.line.set_ydata([event.xdata, pow(10, (np.log10(event.xdata) + np.log10(1000000) - np.log10(0.0001)))])
            self.c.draw_idle()

        def releaseonclick(self, event):
            global k, s, c
            if self.o == "h":
                self.XorY = self.line.get_ydata()[0]
                if params["test_type"] == "BU":
                    m = self.XorY
                    k = 70.6 * params["qo"] * params["bo"] * params["muo"] / (params["h"] * m)
                    print("Estimated Permeability is", round(k, 4), "md")
                    s = 1.151 * ((params["Pi"] - params["pwf"]) / (2.303 * m) - np.log10(k / (params["PHIE"] * params["muo"] * params["ct"] * params["rw"] ** 2)) + 3.23 - np.log10(params["tp"]))
                    print("Estimated skin is", round(s, 2))
                else:
                    m = self.XorY
                    k = 70.6 * params["qo"] * params["bo"] * params["muo"] / (params["h"] * m)
                    print("Estimated Permeability is", round(k, 4), "md")
                    t_fit = derv_data.loc[(derv_data["derv"] >= m*0.9) & (derv_data["derv"] <= m*1.1), "t"].median()
                    pwf_fit = derv_data.loc[(derv_data["derv"] >= m * 0.9) & (derv_data["derv"] <= m * 1.1), "p"].median()
                    s = 1.151 * ((params["Pi"] - pwf_fit) / (2.303 * m) - np.log10(k * t_fit / (params["PHIE"] * params["muo"] * params["ct"] * params["rw"] ** 2)) + 3.23)
                    print("Estimated skin is", round(s, 2))

            elif self.o == "v":
                self.XorY = self.line.get_xdata()[0]
            elif self.o == "US":
                self.XorY = self.line.get_ydata()[0]
                Der = self.XorY
                if params["test_type"] == "BU":
                    c = params["qo"] * params["bo"] * 0.0001 / (24 * Der)
                    print("Wellbore Storage,", round(c, 5))
                else:
                    c = params["qo"] * params["bo"] * 0.0001 / (24 * Der)
                    print("Wellbore Storage,", round(c, 5))

            self.c.mpl_disconnect(self.releaser)
            self.c.mpl_disconnect(self.follower)

    # plot Figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(derv_data["t"], derv_data["dp"], color="blue", linestyle="none", marker="o", label="Delta Pressure")
    plt.plot(derv_data["t"], derv_data["derv"], color="red", linestyle="none", marker="o", label="Pressure Derivative")
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(visible=True, which='both')
    k_line = draggable_lines(ax, "h", np.mean(derv_data["derv"]))
    US_line = draggable_lines(ax, "US", 0.01)
    plt.suptitle("Pressure Analysis Plot")
    plt.title("Move the lines then click to estimate parameters")
    plt.show()
    try:
        return k, s, c
    except NameError:
        try:
            return k, s, np.nan
        except NameError:
            print("--------------------------")
            print("Parameters were not estimated")
            print("--------------------------")
            return np.nan, np.nan, np.nan


# read raw pressure data file
p_data = pd.read_csv("p_data.csv")
p_data["DateTime"] = pd.to_datetime(p_data["Date"] + " " + p_data["Time"])
# Define Input parameters
params_dict = {"bo": 1.5, "muo": 0.35, "qo": 800, "h": 40, "PHIE": 0.12, "Pi": 5410, "ct": 1e-5, "rw": 0.3}


DD_data, BU_data = get_limits(p_data)
data, params_dict = prepare_data(DD_data, BU_data, params_dict, "BU")
final_data, params_dict = calc_der(data, params_dict, 0.5)
k, s, c = derivative_plot_analysis(final_data, params_dict)

print("--------------------------")
print("Draw-down Duration", round(params_dict["tp"], 2), "hr")
print("Draw-down Rate", params_dict["qo"], "STB/D")
print("Buildup Duration", round(params_dict["BU_duration"], 2), "hr")
print("Analyzed test", params_dict["test_type"])
print("Final Permeability is {}, md".format(round(k, 5)))
print("Final Skin is {}".format(round(s, 5)))
print("Final Wellbore Storage is {}, bbl/psi".format(round(c, 5)))
print("--------------------------")
