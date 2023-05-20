const chart_speed = new Chart(document.getElementById("speed").getContext("2d"), {
    data: {
        datasets: [
            {type: "bar", label: "MS1", stack: "1"},
            {type: "bar", label: "MS2", stack: "1"},
            {type: "line", label: `average MS1 (${F1_MEAN.toFixed(2)} Hz)`, pointStyle: false},
            {type: "line", label: `average MS2 (${F2_MEAN.toFixed(2)} Hz)`, pointStyle: false},
        ],
    },
    options: {
        scales: {
            x: {title: {display: true, text: "retention time (min)"}},
            y: {
                title: {display: true, text: "acquisition speed (Hz)"},
                ticks: {callback: (v, i, t) => v + " Hz"},
                beginAtZero: true,
            },
        },
    },
})

function plot_speed(c=chart_speed) {
    var bin = document.getElementById("speed_bin").value
    var n_bin = Math.ceil(RT_MAX / bin) + 1
    var y1 = Array(n_bin).fill(0)
    var y2 = Array(n_bin).fill(0)
    RT1.forEach(x => y1[Math.floor(x / bin)] += 1 / bin)
    RT2.forEach(x => y2[Math.floor(x / bin)] += 1 / bin)
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat((i * bin / 60).toFixed(4)))
    c.data.datasets[0].data = y1
    c.data.datasets[1].data = y2
    c.data.datasets[2].data = Array(n_bin).fill(F1_MEAN)
    c.data.datasets[3].data = Array(n_bin).fill(F2_MEAN)
    c.update()
}

const chart_mz = new Chart(document.getElementById("mz").getContext("2d"), {
    data: {labels: RT2.map(x => x / 60), datasets: [{type: "scatter", label: "MS2", data: MZ2}]},
    options: {
        scales: {
            x: {title: {display: true, text: "retention time (min)"}, min: 0, max: Math.ceil(RT_MAX / 60)},
            y: {title: {display: true, text: "activation center (Th)"}, min: Math.floor(MZ2_MIN / 100) * 100, max: Math.ceil(MZ2_MAX / 100) * 100},
        },
        elements: {point: {radius: 1}},
        animation: false,
        showLine: false,
        spanGaps: true,
        plugins: {legend: {display: false}},
    },
})

const chart_it = new Chart(document.getElementById("it").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "MS1"}, {type: "bar", label: "MS2"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "ion injection time (ms)"}},
            y: {title: {display: true, text: "#scan"}, beginAtZero: true},
        },
    },
})

function plot_it(c=chart_it) {
    var min_rt = document.getElementById("it_min_rt").value * 60
    var max_rt = document.getElementById("it_max_rt").value * 60
    var bin = document.getElementById("it_bin").value
    var n_bin = Math.ceil(IT_MAX / bin) + 1
    var y1 = Array(n_bin).fill(0)
    var y2 = Array(n_bin).fill(0)
    RT1.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y1[Math.floor(IT1[i] / bin)] += 1})
    RT2.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y2[Math.floor(IT2[i] / bin)] += 1})
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat((i * bin).toFixed(4)))
    c.data.datasets[0].data = y1
    c.data.datasets[1].data = y2
    c.update()
}

const chart_tic = new Chart(document.getElementById("tic").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "MS1"}, {type: "bar", label: "MS2"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "log₁₀(total ion current + 1)"}},
            y: {title: {display: true, text: "#scan"}, beginAtZero: true},
        },
    },
})

function plot_tic(c=chart_tic) {
    var min_rt = document.getElementById("tic_min_rt").value * 60
    var max_rt = document.getElementById("tic_max_rt").value * 60
    var bin = document.getElementById("tic_bin").value
    var xmin = Math.floor(TIC_MIN / bin)
    var n_bin = Math.ceil(TIC_MAX / bin) + 1 - xmin
    var y1 = Array(n_bin).fill(0)
    var y2 = Array(n_bin).fill(0)
    RT1.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y1[Math.floor(TIC1[i] / bin) - xmin] += 1})
    RT2.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y2[Math.floor(TIC2[i] / bin) - xmin] += 1})
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat(((i + xmin) * bin).toFixed(4)))
    c.data.datasets[0].data = y1
    c.data.datasets[1].data = y2
    c.update()
}

const chart_bpi = new Chart(document.getElementById("bpi").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "MS1"}, {type: "bar", label: "MS2"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "log₁₀(base peak intensity + 1)"}},
            y: {title: {display: true, text: "#scan"}, beginAtZero: true},
        },
    },
})

function plot_bpi(c=chart_bpi) {
    var min_rt = document.getElementById("bpi_min_rt").value * 60
    var max_rt = document.getElementById("bpi_max_rt").value * 60
    var bin = document.getElementById("bpi_bin").value
    var xmin = Math.floor(BPI_MIN / bin)
    var n_bin = Math.ceil(BPI_MAX / bin) + 1 - xmin
    var y1 = Array(n_bin).fill(0)
    var y2 = Array(n_bin).fill(0)
    RT1.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y1[Math.floor(BPI1[i] / bin) - xmin] += 1})
    RT2.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y2[Math.floor(BPI2[i] / bin) - xmin] += 1})
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat(((i + xmin) * bin).toFixed(4)))
    c.data.datasets[0].data = y1
    c.data.datasets[1].data = y2
    c.update()
}

const chart_bpm = new Chart(document.getElementById("bpm").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "MS1"}, {type: "bar", label: "MS2"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "base peak mass (Th)"}},
            y: {title: {display: true, text: "#scan"}, beginAtZero: true},
        },
    },
})

function plot_bpm(c=chart_bpm) {
    var min_rt = document.getElementById("bpm_min_rt").value * 60
    var max_rt = document.getElementById("bpm_max_rt").value * 60
    var bin = document.getElementById("bpm_bin").value
    var xmin = Math.floor(BPM_MIN / bin)
    var n_bin = Math.ceil(BPM_MAX / bin) + 1 - xmin
    var y1 = Array(n_bin).fill(0)
    var y2 = Array(n_bin).fill(0)
    RT1.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y1[Math.floor(BPM1[i] / bin) - xmin] += 1})
    RT2.forEach((v, i) => {if (v >= min_rt && v <= max_rt) y2[Math.floor(BPM2[i] / bin) - xmin] += 1})
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat(((i + xmin) * bin).toFixed(4)))
    c.data.datasets[0].data = y1
    c.data.datasets[1].data = y2
    c.update()
}

document.getElementById("it_max_rt").value = Math.ceil(RT_MAX / 60)
document.getElementById("tic_max_rt").value = Math.ceil(RT_MAX / 60)
document.getElementById("bpi_max_rt").value = Math.ceil(RT_MAX / 60)
document.getElementById("bpm_max_rt").value = Math.ceil(RT_MAX / 60)

plot_speed()
plot_it()
plot_tic()
plot_bpi()
plot_bpm()
