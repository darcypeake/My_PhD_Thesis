import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.backends.backend_pdf import PdfPages

def project(theta, phi):
    x = theta
    y = np.sin(theta) * (phi - np.pi) / np.pi
    return x, y

plot_paper = True

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
      r"$k_t$":      ["kt2Pi", "ktPi", "kt", "kt08", "ktPio3", "ktPio5"],
      r"C/A":        ["ca2Pi", "caPi", "ca", "ca08", "caPio3", "caPio5"],
      r"anti-$k_t$": ["akt2Pi","aktPi","akt","akt08","aktPio3", "aktPio5"],
      r"Centauro":   ["cen2Pi", "cenPi","cen",  "cen08",  "cenPio3", "cenPio5"],
    }

    titles = {
        "kt2Pi": r"($R_0=4$)",
        "ktPi":  r"($R_0=2$)",
        "kt":    r"($R_0=1$)",
        "kt08":  r"($R_0=0.8$)",
        "ktPio3":  r"($R_0=2/3$)",
        "ktPio5":  r"($R_0=2/5$)",

        "ca2Pi": r"($R_0=4$)",
        "caPi":  r"($R_0=2$)",
        "ca":    r"($R_0=1$)",
        "ca08":  r"($R_0=0.8$)",
        "caPio3":  r"($R_0=2/3$)",
        "caPio5":  r"($R_0=2/5$)",

        "akt2Pi": r"($R_0=4$)",
        "aktPi":  r"($R_0=2$)",
        "akt":    r"($R_0=1$)",
        "akt08":  r"($R_0=0.8$)",
        "aktPio3":  r"($R_0=2/3$)",
        "aktPio5":  r"($R_0=2/5$)",

        "cen2Pi": r"($R_0=4$)",
        "cenPi":  r"($R_0=2$)",
        "cen":    r"($R_0=1$)",
        "cen08":  r"($R_0=0.8$)",
        "cenPio3":  r"($R_0=2/3$)",
        "cenPio5":  r"($R_0=2/5$)"
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
                            alpha=alpha,
                            rasterized=True)

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
      
def plot_event_paper(inputs, pdfname="event.pdf"):
    algo_groups = {
      r"$k_t$":      ["kt2Pi", "ktPi", "kt", "kt08", "ktPio3", "ktPio5"],
      r"C/A":        ["ca2Pi", "caPi", "ca", "ca08", "caPio3", "caPio5"],
      r"anti-$k_t$": ["akt2Pi","aktPi","akt","akt08","aktPio3", "aktPio5"],
      r"Centauro":   ["cen2Pi", "cenPi","cen",  "cen08",  "cenPio3", "cenPio5"],
    }

    titles = {
        "kt2Pi":   r"$R=\frac{\pi}{2}R_0=2\pi$",
        "ktPi":    r"$R=\frac{\pi}{2}R_0=\pi$",
        "kt":      r"$R=\frac{\pi}{2}R_0=\frac{\pi}{2}$",
        "kt08":    r"$R=\frac{\pi}{2}R_0=0.4\pi$",
        "ktPio3":  r"$R=\frac{\pi}{2}R_0=\frac{\pi}{3}$",
        "ktPio5":  r"$R=\frac{\pi}{2}R_0=\frac{\pi}{5}$",

        "ca2Pi":   r"$R=\frac{\pi}{2}R_0=2\pi$",
        "caPi":    r"$R=\frac{\pi}{2}R_0=\pi$",
        "ca":      r"$R=\frac{\pi}{2}R_0=\frac{\pi}{2}$",
        "ca08":    r"$R=\frac{\pi}{2}R_0=0.4\pi$",
        "caPio3":  r"$R=\frac{\pi}{2}R_0=\frac{\pi}{3}$",
        "caPio5":  r"$R=\frac{\pi}{2}R_0=\frac{\pi}{5}$",

        "akt2Pi":   r"$R=\frac{\pi}{2}R_0=2\pi$",
        "aktPi":    r"$R=\frac{\pi}{2}R_0=\pi$",
        "akt":      r"$R=\frac{\pi}{2}R_0=\frac{\pi}{2}$",
        "akt08":    r"$R=\frac{\pi}{2}R_0=0.4\pi$",
        "aktPio3":  r"$R=\frac{\pi}{2}R_0=\frac{\pi}{3}$",
        "aktPio5":  r"$R=\frac{\pi}{2}R_0=\frac{\pi}{5}$",

        "cen2Pi":   r"$R=R_0=4$",
        "cenPi":    r"$R=R_0=2$",
        "cen":      r"$R=R_0=1$",
        "cen08":    r"$R=R_0=0.8$",
        "cenPio3":  r"$R=R_0=\frac{2}{3}$",
        "cenPio5":  r"$R=R_0=\frac{2}{5}$"
    }
    
    njets_for_set = 0
    with PdfPages(pdfname) as pdf: 
      for set in inputs.keys():
        print(set)
        iset = float(str(set).replace("set",""))
        # if (iset % 2 != 0):
        #     continue
        input = inputs[set]
        xdis = input["x"]
        Q    = input["Q"]
        ydis = input["y"]
        # create the figure 
        fig, axs = plt.subplots(1, 4, figsize=(16, 4))
        for j, (name, label) in enumerate(zip("ktPio3 caPio3 aktPio3 cenPio3".split(), ["$k_t$", "C/A", "anti-$k_t$", "Centauro"])):
            
            particles, labels = convert(input[name]["data"])
            #print(name, particles)
            theta  = np.array([p["theta"] for p in particles])
            phi    = np.array([p["phi"]   for p in particles])
            mom    = np.array([p["E"]     for p in particles])
            labels = np.array(labels)

            macro_label = input[name]["mlabel"]

            x, y = project(theta, phi)

            norm = (mom - mom.min()) / (mom.max() - mom.min() + 1e-9)
            sizes = np.where(mom == 1e-6, 5, 20 + 1000 * norm**0.3)

            cmap = plt.get_cmap("terrain")
            ax = axs[j]

            njets = np.max(labels) + 1
            if (iset % 2 != 0):
                njets_for_set = njets
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
                        alpha=alpha,
                        rasterized=True)

            # ---- FORMATTING ----
            ax.set_xlim(-0.2, np.pi+0.2)
            ax.set_ylim(-1.1, 1.1)
            ax.set_axis_off()

            ax.text(0.01, 0.995, 'Remnant', transform=ax.transAxes,
                    ha='left', va='top', fontsize=10)
            ax.text(0.99, 0.985, 'Current', transform=ax.transAxes,
                    ha='right', va='top', fontsize=10)

            ax.text(1.00, 0.05, titles[name], fontsize = 12,
                    transform=ax.transAxes, ha="right")

            ax.text(-.0, 0.05,
                    f"{label}", fontsize = 12,
                    transform=ax.transAxes)

        # fig.suptitle(
        #     rf"$x_{{\rm B}}$={round(xdis,3)}, $Q$={round(Q,1)} GeV",
        #     fontsize=14,
        #     y=0.98
        # )

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

import ast
import glob

if(not plot_paper):
    input_data = {}
    set_counter = 1

    # for filename in sorted(glob.glob("../../panscales/plot_akt*.dat")):
    for filename in sorted(glob.glob("input_data_Rdep/plot_*.dat")):
        with open(filename) as f:
            print(f"Reading {filename}...")
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

    plot_event(input_data, "events_Rdep.pdf")
else:
    input_data = {}
    set_counter = 1
    for filename in sorted(glob.glob("input_data_Rdep/plot_16.dat")):
        with open(filename) as f:
            print(f"Reading {filename}...")
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
    plot_event_paper(input_data, "events_plot_paper.pdf")