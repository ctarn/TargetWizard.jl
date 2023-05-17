const chart_speed = new Chart(document.getElementById("speed").getContext("2d"), {
    data: {
        datasets: [
            {type: "bar", label: "MS1", stack: "1"},
            {type: "bar", label: "MS2", stack: "1"},
            {type: "line", label: `average MS1 (${f1_mean.toFixed(2)} Hz)`, pointStyle: false},
            {type: "line", label: `average MS2 (${f2_mean.toFixed(2)} Hz)`, pointStyle: false},
        ],
    },
    options: {
        scales: {
            x: {title: {display: true, text: "retention time (min)"}},
            y: {
                title: {display: true, text: "acquisition speed (Hz)"},
                ticks: {callback: (v, i, t) => v + " Hz"},
            },
        },
    },
})

function plot_speed() {
    bin = document.getElementById("speed_bin").value
    n_bin = Math.ceil(rt_max / bin) + 1
    f1 = Array(n_bin).fill(0)
    f2 = Array(n_bin).fill(0)
    rt1.forEach(x => f1[Math.floor(x / bin)] += 1 / bin)
    rt2.forEach(x => f2[Math.floor(x / bin)] += 1 / bin)
    chart_speed.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat((i * bin / 60).toFixed(4)))
    chart_speed.data.datasets[0].data = f1
    chart_speed.data.datasets[1].data = f2
    chart_speed.data.datasets[2].data = Array(n_bin).fill(f1_mean)
    chart_speed.data.datasets[3].data = Array(n_bin).fill(f2_mean)
    chart_speed.update()
}

const chart_mz = new Chart(document.getElementById("mz").getContext("2d"), {
    data: {labels: rt2.map(x => x / 60), datasets: [{type: "scatter", label: "MS2", data: mz2}]},
    options: {
        scales: {
            x: {title: {display: true, text: "retention time (min)"}, min: 0, max: Math.ceil(rt_max / 60)},
            y: {title: {display: true, text: "activation center (Th)"}, min: Math.floor(mz2_min / 100) * 100, max: Math.ceil(mz2_max / 100) * 100},
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
            y: {title: {display: true, text: "#scan"}},
        },
    },
})

function plot_it() {
    lower = document.getElementById("it_min_rt").value * 60
    upper = document.getElementById("it_max_rt").value * 60
    bin = document.getElementById("it_bin").value
    n_bin = Math.ceil(it_max / bin) + 1
    f1 = Array(n_bin).fill(0)
    f2 = Array(n_bin).fill(0)
    rt1.forEach((v, i) => {if (v >= lower && v <= upper) f1[Math.floor(it1[i] / bin)] += 1})
    rt2.forEach((v, i) => {if (v >= lower && v <= upper) f2[Math.floor(it2[i] / bin)] += 1})
    chart_it.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat((i * bin).toFixed(4)))
    chart_it.data.datasets[0].data = f1
    chart_it.data.datasets[1].data = f2
    chart_it.update()
}

const chart_tic = new Chart(document.getElementById("tic").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "MS1"}, {type: "bar", label: "MS2"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "log₁₀(total ion current + 1)"}},
            y: {title: {display: true, text: "#scan"}},
        },
    },
})

function plot_tic() {
    lower = document.getElementById("tic_min_rt").value * 60
    upper = document.getElementById("tic_max_rt").value * 60
    bin = document.getElementById("tic_bin").value
    m = Math.floor(tic_min / bin)
    n_bin = Math.ceil(tic_max / bin) + 1 - m
    f1 = Array(n_bin).fill(0)
    f2 = Array(n_bin).fill(0)
    rt1.forEach((v, i) => {if (v >= lower && v <= upper) f1[Math.floor(tic1[i] / bin) - m] += 1})
    rt2.forEach((v, i) => {if (v >= lower && v <= upper) f2[Math.floor(tic2[i] / bin) - m] += 1})
    chart_tic.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat(((i + m) * bin).toFixed(4)))
    chart_tic.data.datasets[0].data = f1
    chart_tic.data.datasets[1].data = f2
    chart_tic.update()
}

const chart_bpi = new Chart(document.getElementById("bpi").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "MS1"}, {type: "bar", label: "MS2"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "log₁₀(base peak intensity + 1)"}},
            y: {title: {display: true, text: "#scan"}},
        },
    },
})

function plot_bpi() {
    lower = document.getElementById("bpi_min_rt").value * 60
    upper = document.getElementById("bpi_max_rt").value * 60
    bin = document.getElementById("bpi_bin").value
    m = Math.floor(bpi_min / bin)
    n_bin = Math.ceil(bpi_max / bin) + 1 - m
    f1 = Array(n_bin).fill(0)
    f2 = Array(n_bin).fill(0)
    rt1.forEach((v, i) => {if (v >= lower && v <= upper) f1[Math.floor(bpi1[i] / bin) - m] += 1})
    rt2.forEach((v, i) => {if (v >= lower && v <= upper) f2[Math.floor(bpi2[i] / bin) - m] += 1})
    chart_bpi.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat(((i + m) * bin).toFixed(4)))
    chart_bpi.data.datasets[0].data = f1
    chart_bpi.data.datasets[1].data = f2
    chart_bpi.update()
}

const chart_bpm = new Chart(document.getElementById("bpm").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "MS1"}, {type: "bar", label: "MS2"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "base peak mass (Th)"}},
            y: {title: {display: true, text: "#scan"}},
        },
    },
})

function plot_bpm() {
    lower = document.getElementById("bpm_min_rt").value * 60
    upper = document.getElementById("bpm_max_rt").value * 60
    bin = document.getElementById("bpm_bin").value
    m = Math.floor(bpm_min / bin)
    n_bin = Math.ceil(bpm_max / bin) + 1 - m
    f1 = Array(n_bin).fill(0)
    f2 = Array(n_bin).fill(0)
    rt1.forEach((v, i) => {if (v >= lower && v <= upper) f1[Math.floor(bpm1[i] / bin) - m] += 1})
    rt2.forEach((v, i) => {if (v >= lower && v <= upper) f2[Math.floor(bpm2[i] / bin) - m] += 1})
    chart_bpm.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat(((i + m) * bin).toFixed(4)))
    chart_bpm.data.datasets[0].data = f1
    chart_bpm.data.datasets[1].data = f2
    chart_bpm.update()
}

document.getElementById("it_max_rt").value = Math.ceil(rt_max / 60)
document.getElementById("tic_max_rt").value = Math.ceil(rt_max / 60)
document.getElementById("bpi_max_rt").value = Math.ceil(rt_max / 60)
document.getElementById("bpm_max_rt").value = Math.ceil(rt_max / 60)

plot_speed()
plot_it()
plot_tic()
plot_bpi()
plot_bpm()
