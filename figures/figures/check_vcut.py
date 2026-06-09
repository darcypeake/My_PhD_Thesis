import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.backends.backend_pdf import PdfPages

def project(theta, phi):
    x = theta
    y = np.sin(theta) * (phi - np.pi) / np.pi
    return x, y


# def plot_event(particles, labels, title="Event display", pdf="event.pdf"):
    
# ------------------ conversion ------------------

def convert(data):
    particles = []
    labels = []
    for E, theta, phi, lab in data:
        particles.append({"E":E, "theta":theta, "phi":phi})
        labels.append(lab)
    return particles, labels

# plot
def plot_event(inputs, pdfname="event.pdf"):
    algo_groups = {
      r"$k_t$":  ["ktycut1", "ktycut0.1", "ktycut0.01", "ktycut0.001", "ktycut0.0001"],
      r"C/A":  ["caycut1", "caycut0.1", "caycut0.01", "caycut0.001", "caycut0.0001"],
    }
          
    titles = {
        "ktycut1": r"($y_{\rm cut}=1$)",
        "ktycut0.1": r"($y_{\rm cut}=0.1$)",
        "ktycut0.01": r"($y_{\rm cut}=0.01$)",
        "ktycut0.001": r"($y_{\rm cut}=0.001$)",
        "ktycut0.0001": r"($y_{\rm cut}=0.0001$)",

        "caycut1": r"($y_{\rm cut}=1$)",
        "caycut0.1": r"($y_{\rm cut}=0.1$)",
        "caycut0.01": r"($y_{\rm cut}=0.01$)",
        "caycut0.001": r"($y_{\rm cut}=0.001$)",
        "caycut0.0001": r"($y_{\rm cut}=0.0001$)",
    }
    
    njets_for_set = 0
    with PdfPages(pdfname) as pdf: 
      for set in inputs.keys():
        print(set)
        iset = float(str(set).replace("set",""))
        input = inputs[set]
        xdis = input["x"]
        Q    = input["Q"]
        ydis = input["y"]
        # create the figure 
        for algo, group_names in algo_groups.items():
          fig, axs = plt.subplots(1, len(group_names),
                                  figsize=(16/4*len(group_names), 4))

          njets = 0
          for j, name in enumerate(group_names):
              particles, labels = convert(input[name]["data"])
              #print(name, particles)
              theta  = np.array([p["theta"] for p in particles])
              phi    = np.array([p["phi"]   for p in particles])
              mom    = np.array([p["E"]     for p in particles])
              labels = np.array(labels)

              macro_label = input[name]["mlabel"]

              njets = np.max(labels) + 1
              if (iset % 2 != 0):
                  njets_for_set = njets

              x, y = project(theta, phi)

              norm = (mom - mom.min()) / (mom.max() - mom.min() + 1e-9)
              sizes = np.where(mom == 1e-6, 5, 20 + 1000 * norm**0.3)

              cmap = plt.get_cmap("terrain")
              ax = axs[j]

              # ---- GRID ----
              for t in np.linspace(0, np.pi, 7):
                  ph = np.linspace(0, 2*np.pi, 300)
                  Xg, Yg = project(np.full_like(ph, t), ph)
                  ax.plot(Xg, Yg, color='0.85', lw=1)

              for ph0 in [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]:
                  th = np.linspace(0, np.pi, 300)
                  Xg, Yg = project(th, np.full_like(th, ph0))
                  ax.plot(Xg, Yg, color='0.85', lw=1)

              # ---- PARTICLES ----
              ncol_lab = 0
              for lab in np.unique(labels):
                  m = labels == lab
                  all_soft = np.all(np.isclose(mom[m], 1e-6))

                  if lab == macro_label:
                      c = 'red'
                      ec = 'black'
                      lw = 0.1
                      alpha = 0.85
                  else:
                      if all_soft:
                          normc = lab / max(1, njets - 1)
                          c = plt.cm.Greys(normc)
                          ec = 'gray'
                          lw = 0.3
                          alpha = 0.5
                      else:
                          normc = ncol_lab / max(1, njets_for_set - 1)
                          c = cmap(normc)
                          ec = 'gray'
                          lw = 0.2
                          alpha = 0.6
                          ncol_lab += 1

                  ax.scatter(x[m], y[m],
                            s=sizes[m],
                            color=c,
                            edgecolors=ec,
                            linewidth=lw,
                            alpha=alpha)

              # ---- FORMATTING ----
              ax.set_xlim(-0.2, np.pi+0.2)
              ax.set_ylim(-1.1, 1.1)
              ax.set_axis_off()

              ax.text(0.01, 0.995, 'Remnant', transform=ax.transAxes,
                      ha='left', va='top', fontsize=8)
              ax.text(0.99, 0.985, 'Current', transform=ax.transAxes,
                      ha='right', va='top', fontsize=8)

              ax.text(-.05, 0.05, titles[name],
                      transform=ax.transAxes)

              ax.text(1.05, 0.05,
                      f"$n_{{\\rm jets}}$ ={njets}",
                      transform=ax.transAxes, ha="right")

          fig.suptitle(
              rf"{algo} jets, $x_{{\rm B}}$={round(xdis,3)}, $Q$={round(Q,1)} GeV",
              fontsize=14,
              y=0.98
          )

          pdf.savefig(fig, bbox_inches='tight')
          plt.close(fig)
      
        

import ast
import glob

input_data = {}
set_counter = 1

# for filename in sorted(glob.glob("../../panscales/plot_akt*.dat")):
for filename in sorted(glob.glob("input_data_vcutdep/plot_*.dat")):
    with open(filename) as f:
        # Remove comments and empty lines
        lines = [
            line.split("#", 1)[0]
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]

    text = "".join(lines)
    text = "{" + text + "}"  # make it a full dict

    data = ast.literal_eval(text)

    # Loop through sets in this file
    for _, value in data.items():
        new_key = f"set{set_counter}"
        input_data[new_key] = value
        set_counter += 1

# with open("../../panscales/plot.dat") as f:
#   text = f.read()

# # Ensure it's a valid Python dict
# text = "{" + text + "}"

# input_data = ast.literal_eval(text)
plot_event(input_data, "events_vcutdep.pdf")
#################
## xdis = 0.468, Q = 87.946, y = 0.183
#################
# kt_data = [
# (2.327, 2.826, 2.102, 0),(16.974, 3.089, 5.057, 0),(3.443, 3.095, 5.276, 0),(1.992, 2.855, 5.005, 0),(15.715, 3.057, 5.638, 0),(2.553, 2.895, 6.063, 0),(0.636, 1.475, 0.661, 0),(0.836, 1.737, 6.054, 0),(0.527, 1.074, 1.994, 1),(0.804, 0.923, 2.725, 1),(0.556, 0.636, 2.220, 1),(1.565, 0.394, 1.836, 1),(0.927, 0.464, 1.026, 1),(2.231, 0.268, 1.770, 1),(10.384, 0.055, 3.340, 2),(25.681, 0.041, 3.848, 2),(6.740, 0.053, 1.820, 3)
# ]
# ca_data = [
# (2.327, 2.826, 2.102, 0),(1.992, 2.855, 5.005, 0),(2.553, 2.895, 6.063, 0),(15.715, 3.057, 5.638, 0),(16.974, 3.089, 5.057, 0),(3.443, 3.095, 5.276, 0),(0.636, 1.475, 0.661, 0),(0.836, 1.737, 6.054, 0),(0.527, 1.074, 1.994, 1),(0.556, 0.636, 2.220, 1),(0.804, 0.923, 2.725, 1),(0.927, 0.464, 1.026, 2),(1.565, 0.394, 1.836, 2),(2.231, 0.268, 1.770, 2),(6.740, 0.053, 1.820, 3),(10.384, 0.055, 3.340, 4),(25.681, 0.041, 3.848, 4)
# ]
# centauro_data = [
# (1.565, 0.394, 1.836, 0),(6.740, 0.053, 1.820, 0),(25.681, 0.041, 3.848, 1),(0.927, 0.464, 1.026, 2),(10.384, 0.055, 3.340, 3),(2.231, 0.268, 1.770, 4),(0.636, 1.475, 0.661, 5),(0.804, 0.923, 2.725, 6),(0.556, 0.636, 2.220, 7),(0.527, 1.074, 1.994, 7),(0.836, 1.737, 6.054, 8),(2.327, 2.826, 2.102, 8),(2.553, 2.895, 6.063, 8),(1.992, 2.855, 5.005, 8),(15.715, 3.057, 5.638, 8),(16.974, 3.089, 5.057, 8),(3.443, 3.095, 5.276, 8)
# ]

################
## xdis = 0.257, Q = 44.669, y = 0.086
################
# kt_data = [
# (1.055, 2.392, 0.429, 0),(2.621, 2.808, 0.684, 0),(0.752, 2.888, 1.643, 0),(1.935, 2.721, 4.218, 0),(2.195, 2.691, 3.378, 0),(0.393, 2.233, 1.848, 0),(1.468, 2.717, 2.598, 0),(0.492, 2.757, 2.160, 0),(5.502, 3.040, 5.582, 0),(5.771, 3.064, 4.181, 0),(1.173, 0.861, 5.454, 1),(0.569, 1.635, 5.607, 1),(0.203, 1.333, 6.014, 1),(3.419, 0.187, 1.820, 2),(27.317, 0.022, 3.875, 3),(6.355, 0.092, 0.235, 4),(2.860, 0.079, 1.812, 5),(22.722, 0.012, 2.200, 5)
# ]
# ca_data = [
# (0.393, 2.233, 1.848, 0),(1.055, 2.392, 0.429, 0),(1.468, 2.717, 2.598, 0),(0.492, 2.757, 2.160, 0),(5.502, 3.040, 5.582, 0),(5.771, 3.064, 4.181, 0),(2.621, 2.808, 0.684, 0),(0.752, 2.888, 1.643, 0),(1.935, 2.721, 4.218, 0),(2.195, 2.691, 3.378, 0),(1.173, 0.861, 5.454, 1),(0.569, 1.635, 5.607, 1),(0.203, 1.333, 6.014, 1),(3.419, 0.187, 1.820, 2),(6.355, 0.092, 0.235, 3),(2.860, 0.079, 1.812, 4),(27.317, 0.022, 3.875, 5),(22.722, 0.012, 2.200, 6)
# ]
# centauro_data = [
# (2.860, 0.079, 1.812, 0),(3.419, 0.187, 1.820, 0),(0.393, 2.233, 1.848, 1),(22.722, 0.012, 2.200, 2),(6.355, 0.092, 0.235, 3),(27.317, 0.022, 3.875, 4),(0.203, 1.333, 6.014, 5),(0.569, 1.635, 5.607, 5),(1.173, 0.861, 5.454, 5),(1.055, 2.392, 0.429, 6),(1.468, 2.717, 2.598, 6),(0.492, 2.757, 2.160, 6),(5.502, 3.040, 5.582, 6),(5.771, 3.064, 4.181, 6),(2.621, 2.808, 0.684, 6),(0.752, 2.888, 1.643, 6),(1.935, 2.721, 4.218, 6),(2.195, 2.691, 3.378, 6)
# ]



##############
## xdis = 0.001, Q = 5.779, y = 0.397
###############
# kt_data = [
# (0.778, 1.754, 3.361, 0),(2.719, 1.632, 3.586, 0),(0.629, 0.896, 0.139, 1),(2.389, 0.165, 1.511, 1),(19.715, 0.034, 0.296, 1),(0.627, 0.213, 6.220, 1),(3.605, 0.103, 0.129, 1),(2.875, 0.083, 5.261, 1),(17.299, 0.051, 5.629, 1),(2.631, 0.081, 6.023, 1),(36.887, 0.007, 6.097, 1),(17.119, 0.009, 0.106, 1),(117.748, 0.002, 0.147, 1),(1481.655, 0.000, 6.196, 1),(365.478, 0.001, 0.201, 1),(314.731, 0.001, 0.312, 1),(3.642, 0.574, 1.868, 2),(0.560, 0.417, 2.527, 2),(1.073, 0.435, 3.798, 3),(10.470, 0.024, 4.030, 3),(138.500, 0.002, 3.418, 3),(492.773, 0.001, 3.779, 3),(58.403, 0.002, 2.114, 4)
# ]
# ca_data = [
# (0.778, 1.754, 3.361, 0),(2.719, 1.632, 3.586, 0),(0.629, 0.896, 0.139, 1),(3.642, 0.574, 1.868, 2),(0.560, 0.417, 2.527, 2),(1.073, 0.435, 3.798, 3),(0.627, 0.213, 6.220, 4),(2.389, 0.165, 1.511, 5),(3.605, 0.103, 0.129, 6),(2.875, 0.083, 5.261, 7),(17.299, 0.051, 5.629, 7),(2.631, 0.081, 6.023, 7),(19.715, 0.034, 0.296, 8),(10.470, 0.024, 4.030, 9),(36.887, 0.007, 6.097, 10),(17.119, 0.009, 0.106, 10),(58.403, 0.002, 2.114, 11),(117.748, 0.002, 0.147, 12),(138.500, 0.002, 3.418, 13),(492.773, 0.001, 3.779, 14),(1481.655, 0.000, 6.196, 15),(365.478, 0.001, 0.201, 15),(314.731, 0.001, 0.312, 15)
# ]
# centauro_data = [
# (2.389, 0.165, 1.511, 0),(0.627, 0.213, 6.220, 1),(138.500, 0.002, 3.418, 2),(36.887, 0.007, 6.097, 3),(10.470, 0.024, 4.030, 4),(58.403, 0.002, 2.114, 5),(365.478, 0.001, 0.201, 6),(117.748, 0.002, 0.147, 7),(314.731, 0.001, 0.312, 8),(492.773, 0.001, 3.779, 9),(1481.655, 0.000, 6.196, 10),(19.715, 0.034, 0.296, 11),(17.299, 0.051, 5.629, 12),(2.875, 0.083, 5.261, 13),(17.119, 0.009, 0.106, 14),(0.560, 0.417, 2.527, 15),(2.631, 0.081, 6.023, 16),(3.605, 0.103, 0.129, 17),(1.073, 0.435, 3.798, 18),(0.778, 1.754, 3.361, 19),(2.719, 1.632, 3.586, 19),(0.629, 0.896, 0.139, 20),(3.642, 0.574, 1.868, 21)
# ]

##############
# xdis = 0.539, Q = 64.307, y = 0.085
##############
# kt_data = [
# (0.505, 1.768, 4.688, 0),(1.084, 1.976, 5.696, 0),(14.140, 3.096, 0.158, 0),(6.936, 3.092, 5.246, 0),(1.235, 2.755, 3.536, 0),(1.418, 2.991, 2.219, 0),(6.819, 3.066, 1.586, 0),(1.033, 3.023, 1.294, 0),(15.480, 0.045, 1.344, 1),(8.777, 0.052, 2.327, 1),(2.281, 0.368, 3.944, 2)
# ]
# ca_data = [
# (1.235, 2.755, 3.536, 0),(1.418, 2.991, 2.219, 0),(6.819, 3.066, 1.586, 0),(1.033, 3.023, 1.294, 0),(14.140, 3.096, 0.158, 0),(6.936, 3.092, 5.246, 0),(0.505, 1.768, 4.688, 0),(1.084, 1.976, 5.696, 0),(2.281, 0.368, 3.944, 1),(8.777, 0.052, 2.327, 2),(15.480, 0.045, 1.344, 3)
# ]
# centauro_data = [
# (1.084, 1.976, 5.696, 0),(2.281, 0.368, 3.944, 1),(15.480, 0.045, 1.344, 2),(0.505, 1.768, 4.688, 3),(8.777, 0.052, 2.327, 4),(1.235, 2.755, 3.536, 5),(14.140, 3.096, 0.158, 5),(6.936, 3.092, 5.246, 5),(1.418, 2.991, 2.219, 5),(6.819, 3.066, 1.586, 5),(1.033, 3.023, 1.294, 5)
# ]


# kt_particles, kt_labels     = convert(kt_data)
# ca_particles, ca_labels     = convert(ca_data)
# cent_particles, cent_labels = convert(centauro_data)

# # ------------------ plotting ------------------

# plot_event(kt_particles, kt_labels, "kT ordering", "kt.pdf")
# plot_event(ca_particles, ca_labels, "CA ordering", "ca.pdf")
# plot_event(cent_particles, cent_labels, "Centauro ordering", "centauro.pdf")