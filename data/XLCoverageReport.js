const chart_{{ id }} = new Chart(document.getElementById("{{ id }}").getContext("2d"), {
    data: {datasets: [{type: "bar", label: "Peptide Pair"}, {type: "bar", label: "Peptide α"}, {type: "bar", label: "Peptide β"}]},
    options: {
        scales: {
            x: {title: {display: true, text: "fragment ion coverage (%)"}},
            y: {title: {display: true, text: "#PSM"}, beginAtZero: true},
        },
    },
})

function plot_{{ id }}(c=chart_{{ id }}) {
    var fdr = document.getElementById("{{ id }}_fdr").value
    var bin = document.getElementById("{{ id }}_bin").value
    var n_bin = Math.ceil(100 / bin) + 1
    var y1 = Array(n_bin).fill(0)
    var y2 = Array(n_bin).fill(0)
    var y3 = Array(n_bin).fill(0)
    FDR.forEach((v, i) => {if (v <= fdr) y1[Math.floor(COV_{{ id }}[i] / bin)] += 1})
    FDR.forEach((v, i) => {if (v <= fdr) y2[Math.floor(COV_A_{{ id }}[i] / bin)] += 1})
    FDR.forEach((v, i) => {if (v <= fdr) y3[Math.floor(COV_B_{{ id }}[i] / bin)] += 1})
    c.data.labels = Array(n_bin).fill(0).map((_, i) => parseFloat((i * bin).toFixed(4)))
    c.data.datasets[0].data = y1
    c.data.datasets[1].data = y2
    c.data.datasets[2].data = y3
    c.update()
}

plot_{{ id }}()
