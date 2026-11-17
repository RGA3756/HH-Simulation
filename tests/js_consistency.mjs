// Headless check of the explorer's JavaScript model (run by export_web.js_consistency).
// Usage: node tests/js_consistency.mjs web/explorer.html 6.3 18.5
// Extracts the <script id="hh-model"> block and the embedded data from the
// shipped HTML file, runs the full membrane experiment at each temperature and
// prints the measurements as JSON.
import fs from "node:fs";
import vm from "node:vm";

const [htmlPath, ...temps] = process.argv.slice(2);
const html = fs.readFileSync(htmlPath, "utf8");
const model = html.match(/<script id="hh-model">([\s\S]*?)<\/script>/)[1];
const data = JSON.parse(html.match(/<script id="hh-data" type="application\/json">([\s\S]*?)<\/script>/)[1]);

const context = vm.createContext({});
vm.runInContext(model, context);
const HH = context.HHModel(data.params);

const out = temps.map(Number).map((T) => {
  const res = HH.run(T, data.amp_primary);
  return { T, r: res.m.r, Q_Na: res.m.Q_Na, Q_min: res.m.Q_min, V_peak: res.m.V_peak,
           spiked: res.spiked, t_end: res.sim.t[res.sim.t.length - 1], V_rest: HH.rest.V };
});
console.log(JSON.stringify(out));
