
document.getElementById("calculate").addEventListener("click",()=>{
  let tm = parseFloat(document.getElementById("tm").value);
  let fnr = parseFloat(document.getElementById("fnr").value);
  let cd = parseFloat(document.getElementById("cd").value);
  let r = parseFloat(document.getElementById("r").value);
  let p = parseFloat(document.getElementById("p").value);
  let tf = parseFloat(document.getElementById("tf").value);
  let n = parseInt(document.getElementById("n").value);
  
  let A = [
    [(r - fnr) / tm, cd],
    [fnr / tm, -cd]
  ];
  let b = [
    [1],
    [fnr / (cd * tm)]
  ];

  // Solução exata
  EIG = numeric.eig(A);
  vals = EIG.lambda.x;
  vects = EIG.E.x;
  let c = numeric.solve(vects, b);
  let h = tf / n;
  let TEMP = numeric.rep([1, n], 0);
  let ESOL = numeric.rep([2, n], 0);
  for(let t = 0; t <= n; t++){
    let psum = 0;
    let csum = 0;
    for(let i = 0; i < vals.length; i++){
      psum = psum + c[i] * vects[0][i] * numeric.exp(vals[i] * h * t);
      csum = csum + c[i] * vects[1][i] * numeric.exp(vals[i] * h * t);
    }
    TEMP[0][t] = h*t;
    ESOL[0][t] = psum * p;
    ESOL[1][t] = csum;
  }
  
  let results = document.getElementById("results");
  results.innerHTML = "";
  let pt = tf/10, ph = 0;
  for(let t = 0; t <= n; t++){
    if (ph <= TEMP[0][t]) {
      ph = ph + pt;
      let tr = document.createElement("tr");
      let td1 = document.createElement("td");
      td1.innerHTML = TEMP[0][t].toFixed(2);
      tr.append(td1);
      let td2 = document.createElement("td");
      td2.innerHTML = ESOL[0][t].toFixed(4);
      tr.append(td2);
      let td3 = document.createElement("td");
      td3.innerHTML = ESOL[1][t].toFixed(4);
      tr.append(td3);
      results.append(tr);
    }
  }

  // Solução numérica
  let met = document.getElementById("met");
  let NSOL = numeric.rep([2, n], 0);
  NSOL[0][0] = b[0][0]; NSOL[1][0] = b[1][0];

  // Método de Euler
  if (met.selectedIndex == 0){
    let k1, colk;
    for(let k = 1; k <= n; k++){
      colk = NSOL.map(x=> x[k-1]);
      k1 = numeric.mul(h, numeric.dotMV(A, colk));
      NSOL[0][k] = NSOL[0][k-1] + k1[0];
      NSOL[1][k] = NSOL[1][k-1] + k1[1];
    }
  }

  // Método de Runge Kutta de ordem 2
  else if (met.selectedIndex == 1){
    let k1, k2, colk;
    for(let k = 1; k <= n; k++){
      colk = NSOL.map(x=> x[k-1]);
      k1 = numeric.mul(h, numeric.dotMV(A, colk));
      k2 = numeric.mul(h, numeric.dotMV(A, numeric.addVV(colk, k1)));
      NSOL[0][k] = NSOL[0][k-1] + (k1[0] + k2[0]) / 2;
      NSOL[1][k] = NSOL[1][k-1] + (k1[1] + k2[1]) / 2;
    }
  }

  // Método de Runge Kutta de ordem 4
  else {
    let k1, k2, k3, k4, colk;
    for(let k = 1; k <= n; k++){
      colk = NSOL.map(x=> x[k-1]);
      k1 = numeric.mul(h, numeric.dotMV(A, colk));
      k2 = numeric.mul(h, numeric.dotMV(A, numeric.addVV(colk, numeric.mul(k1, 0.5))));
      k3 = numeric.mul(h, numeric.dotMV(A, numeric.addVV(colk, numeric.mul(k2, 0.5))));
      k4 = numeric.mul(h, numeric.dotMV(A, numeric.addVV(colk, k3)));
      NSOL[0][k] = NSOL[0][k-1] + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
      NSOL[1][k] = NSOL[1][k-1] + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
    }
  }
  for (let i = 0; i <= n; i++){
    NSOL[0][i] = NSOL[0][i] * p;
  }

  let nresults = document.getElementById("nresults");
  nresults.innerHTML = "";
  ph = 0;
  for(let t = 0; t <= n; t++){
    if (ph <= TEMP[0][t]) {
      ph = ph + pt;
      let tr = document.createElement("tr");
      let td1 = document.createElement("td");
      td1.innerHTML = TEMP[0][t].toFixed(2);
      tr.append(td1);
      let td2 = document.createElement("td");
      td2.innerHTML = NSOL[0][t].toFixed(4);
      tr.append(td2);
      let td3 = document.createElement("td");
      td3.innerHTML = NSOL[1][t].toFixed(4);
      tr.append(td3);
      nresults.append(tr);
    }
  }

  // chart
  document.getElementById("chart").innerHTML = "";
  let trace1 = {
    x: TEMP[0],
    y: ESOL[0],
    mode: 'lines',
    name: "Solução Analítica",
    line: {
      color: "#343a40",
      width: 3,
      dash: "solid"
    }
  };
  let trace2 = {
    x: TEMP[0],
    y: NSOL[0],
    mode: 'lines',
    name: "Solução Numérica",
    line: {
      color: "#dc3545",
      width: 3,
      dash: "dot"
    }
  };
  let data = [ trace1, trace2 ];
  let layout = {
    xaxis: {
      title: "Tempo (s)",
    },
    yaxis: {
      title: "Potência (MW)",
      tickmode: "auto",
      exponentformat: "E",
      showexponent: "all",
      rangemode: "tozero"
    }
  };
  Plotly.newPlot('chart', data, layout);

});

function default_monitor(){
  let monitor = document.querySelector(".monitor");
  monitor.innerHTML = "<img style=' width: 80%; height: 80%; object-fit: contain;' src='monitor.png'>";
}
default_monitor();




