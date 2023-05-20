const chart_z_rt = new Chart(document.getElementById("z_rt").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "z=..."}]},
    options: {
        scales: {
            x: {title: {display: true, text: "retention time (min)"}},
            y: {title: {display: true, text: "#target"}, beginAtZero: true},
        },
    },
})

function plot_z_rt(c=chart_z_rt) {
    var bin = document.getElementById("z_rt_bin").value
    var plot_type = document.getElementById("z_rt_type").value
    var n_bin = Math.ceil(RT_MAX / bin) + 1
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat((i * bin / 60).toFixed(4)))
    c.data.datasets = []
    for (var z = Z_MIN; z <= Z_MAX; z++) {
        var y = Array(n_bin).fill(0)
        if (plot_type == "range") {
            Z.forEach((v, i) => {
                if (v == z)
                    for (j = Math.floor(RT_START[i] / bin); j <= Math.ceil(RT_STOP[i] / bin); j++)
                        y[j] += 1
            })
        } else if (plot_type == "center") {
            Z.forEach((v, i) => {
                if (v == z)
                    y[Math.floor((RT[i]) / bin)] += 1
            })
        }
        c.data.datasets.push({type: "bar", data: y, label: `z=${z}`, stack: "1"})
    }
    c.update()
}

const chart_z_hist = new Chart(document.getElementById("z_hist").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "#target"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "charge state"}},
            y: {title: {display: true, text: "#target"}, beginAtZero: true},
        },
    },
})

function plot_z_hist(c=chart_z_hist) {
    var min_rt = document.getElementById("z_hist_min_rt").value * 60
    var max_rt = document.getElementById("z_hist_max_rt").value * 60
    var n_bin = Z_MAX + 1 - Z_MIN
    c.data.labels = Array(n_bin).fill(0).map((_, i) => i + Z_MIN)
    var y = Array(n_bin).fill(0)
    Z.forEach((v, i) => {if (RT_START[i] >= min_rt && RT_STOP[i] <= max_rt) y[v - Z_MIN] += 1})
    c.data.datasets[0].data = y
    c.data.datasets[0].label = `total: ${y.reduce((a, b) => a + b, 0)}`
    c.update()
}

const chart_mass_rt = new Chart(document.getElementById("mass_rt").getContext("2d"), {
    data: {datasets: [{type: "scatter", label: "z=..."}]},
    options: {
        scales: {
            x: {title: {display: true, text: "retention time (min)"}},
            y: {title: {display: true, text: "mass"}},
        },
        elements: {point: {radius: 2}},
    },
})

function plot_mass_rt(c=chart_mass_rt) {
    var plot_type = document.getElementById("mass_rt_type").value
    c.data.datasets = []
    if (plot_type == "da") {
        var vs = M
        c.options.scales.y.title.text = "mass (Da)"
    } else if (plot_type == "th") {
        var vs = MZ
        c.options.scales.y.title.text = "m/z (Th)"
    }
    for (var z = Z_MIN; z <= Z_MAX; z++) {
        var y = []
        for (var i = 0; i < RT.length; i++) {
            if (Z[i] == z) y.push({x: RT[i] / 60, y: vs[i]})
        }
        c.data.datasets.push({type: "scatter", data: y, label: `z=${z}: ${y.length}`})
    }
    c.update()
}

const chart_mass_hist = new Chart(document.getElementById("mass_hist").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "z=..."}]},
    options: {
        scales: {
            x: {title: {display: true, text: "mass"}},
            y: {title: {display: true, text: "#target"}, beginAtZero: true},
        },
    },
})

function plot_mass_hist(c=chart_mass_hist) {
    var plot_type = document.getElementById("mass_hist_type").value
    var min_rt = document.getElementById("mass_hist_min_rt").value * 60
    var max_rt = document.getElementById("mass_hist_max_rt").value * 60
    var bin = document.getElementById("mass_hist_bin").value
    if (plot_type == "da") {
        var xmin = Math.floor(M_MIN / bin)
        var n_bin = Math.ceil(M_MAX / bin) + 1 - xmin
        var vs = M
        c.options.scales.x.title.text = "mass (Da)"
        document.getElementById("mass_unit").innerText = "Da"
    } else if (plot_type == "th") {
        var xmin = Math.floor(MZ_MIN / bin)
        var n_bin = Math.ceil(MZ_MAX / bin) + 1 - xmin
        var vs = MZ
        c.options.scales.x.title.text = "m/z (Th)"
        document.getElementById("mass_unit").innerText = "Th"
    }
    c.data.datasets = []
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat(((i + xmin) * bin).toFixed(4)))
    for (var z = Z_MIN; z <= Z_MAX; z++) {
        var y = Array(n_bin).fill(0)
        Z.forEach((v, i) => {
            if (v == z && RT_START[i] >= min_rt && RT_STOP[i] <= max_rt)
                y[Math.floor(vs[i] / bin) - xmin] += 1
        })
        c.data.datasets.push({type: "bar", data: y, label: `z=${z}: ${y.reduce((a, b) => a + b, 0)}`, stack: "1"})
    }
    c.update()
}

document.getElementById("z_hist_max_rt").value = Math.ceil(RT_STOP_MAX / 60)
document.getElementById("mass_hist_max_rt").value = Math.ceil(RT_STOP_MAX / 60)

plot_z_rt()
plot_z_hist()
plot_mass_rt()
plot_mass_hist()
