
function onMessage(jsonData) {
    var result = JSON.parse(jsonData);
    document.getElementById("solution").innerText = result.solution;
    document.getElementById("time").innerText = result.time;
}

function solve() {
    toggleProgress();
    var Http = new XMLHttpRequest();
    var url='http://localhost:8080/Gauss_SoLAE_Web_App_war_exploded/gauss/solve';
    Http.open("POST", url, true);
    Http.setRequestHeader("Content-Type", "application/json");
    Http.onreadystatechange= function () {
        if (this.readyState===4 && this.status===200) {
            onMessage(Http.responseText);
        }
        toggleProgress();
    }
    var data = JSON.stringify({
        "matrix": parseMatrix(document.getElementById("matrix").value),
        "mpi": document.getElementById("mpi").value
    });
    Http.send(data);
}

function parseMatrix(data) {
    var lines = data.split("\n");
    var res = new Array(lines.length);
    for (var i=0; i < res.length; i++) {
        var words = lines[i].split(" ");
        res[i] = new Array(words.length);
        for (var j = 0; j < words.length; j++)
            res[i][j] = parseFloat(words[j]);
    }
    return res;
}

function reset() {
    document.getElementById("function").value = "";

    document.getElementById("x1").value = "";
    document.getElementById("x2").value = "";
    document.getElementById("y1").value = "";
    document.getElementById("y2").value = "";
    document.getElementById("z1").value = "";
    document.getElementById("z2").value = "";

    document.getElementById("mpi").value = "";

    document.getElementById("solution").innerText = " ";
    document.getElementById("time").innerText = " ";
}

function toggleProgress() {
    var x = document.getElementById("progress");
    if (x.style.display === "none") {
        x.style.display = "block";
    } else {
        x.style.display = "none";
    }
}
