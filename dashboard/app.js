"use strict";

const benchmarkList = document.getElementById("benchmark-list");
const catalogStatus = document.getElementById("catalog-status");
const selectedModel = document.getElementById("selected-model");
const runButton = document.getElementById("run-button");
const runMessage = document.getElementById("run-message");
const resultsContent = document.getElementById("results-content");
const convergenceBadge = document.getElementById("convergence-badge");
const state = { benchmarks: new Map(), currentId: null, requestId: 0 };

function element(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

function formatEnergy(value) {
  return typeof value === "number" && Number.isFinite(value)
    ? value.toFixed(12)
    : "Unavailable";
}

function formatError(value) {
  return typeof value === "number" && Number.isFinite(value)
    ? value.toExponential(3)
    : "Unavailable";
}

function formatMaybe(value) {
  return typeof value === "number" && Number.isFinite(value)
    ? value.toString()
    : "Unavailable";
}

function resetResults(message, isError = false) {
  resultsContent.hidden = true;
  runMessage.hidden = false;
  runMessage.textContent = message;
  runMessage.classList.toggle("is-error", isError);
  convergenceBadge.className = isError ? "badge bg-red-lt" : "badge bg-secondary-lt";
  convergenceBadge.textContent = isError ? "Error" : "Not run";
}

function makeBenchmarkCard(benchmark) {
  const column = element("div", "col-md-6");
  const card = element("article", "card benchmark-card");
  card.dataset.benchmarkId = benchmark.id;
  const body = element("div", "card-body");
  const top = element("div", "d-flex align-items-center justify-content-between gap-2 mb-2");
  top.append(element("span", "page-pretitle", `MODEL ${benchmark.id === "h2" ? "01" : "02"}`));
  top.append(element("span", "badge bg-secondary-lt", `${benchmark.qubits} qubits`));
  body.append(top);

  const title = element("h3", "card-title", benchmark.formula);
  body.append(title);
  body.append(element("p", "text-secondary mb-2", benchmark.display_name));

  const details = element("dl", "model-meta");
  const rows = [
    ["Geometry", `${benchmark.geometry} (${benchmark.geometry_units})`],
    ["Basis / charge", `${benchmark.basis} / ${benchmark.charge > 0 ? "+" : ""}${benchmark.charge}`],
    ["Electrons / spin", `${benchmark.electron_count} / S = ${benchmark.total_spin}`],
    ["Mapping", benchmark.mapping]
  ];
  for (const [label, value] of rows) {
    const row = element("div");
    row.append(element("dt", "small", label));
    row.append(element("dd", "small", value));
    details.append(row);
  }
  body.append(details);

  const button = element("button", "btn btn-outline-primary w-100", "Select benchmark");
  button.type = "button";
  button.setAttribute("aria-pressed", "false");
  button.addEventListener("click", () => selectBenchmark(benchmark.id));
  body.append(button);
  card.append(body);
  column.append(card);
  return column;
}

function selectBenchmark(id) {
  const benchmark = state.benchmarks.get(id);
  if (!benchmark) return;
  state.requestId += 1;
  state.currentId = id;
  selectedModel.textContent = `${benchmark.display_name} · ${benchmark.basis} · ${benchmark.geometry_units}`;
  runButton.disabled = false;
  runButton.textContent = "Run benchmark";
  resetResults("Ready. Run the selected fixed benchmark to calculate current results.");
  for (const card of benchmarkList.querySelectorAll(".benchmark-card")) {
    const selected = card.dataset.benchmarkId === id;
    card.classList.toggle("is-selected", selected);
    const button = card.querySelector("button");
    button.setAttribute("aria-pressed", selected ? "true" : "false");
    button.textContent = selected ? "Selected" : "Select benchmark";
  }
}

function isFiniteNumber(value) {
  return typeof value === "number" && Number.isFinite(value);
}

function validateResult(result, expectedId) {
  if (!result || !result.benchmark || result.benchmark.id !== expectedId) {
    throw new Error("The benchmark returned data for an unexpected model.");
  }
  const requiredNumericFields = [
    "electronic_energy_hartree",
    "nuclear_repulsion_hartree",
    "total_energy_hartree"
  ];
  const optionalNumericFields = [
    "exact_electronic_energy_hartree",
    "reference_electronic_energy_hartree",
    "variational_error_hartree"
  ];
  if (!requiredNumericFields.every((field) => isFiniteNumber(result[field])) ||
      !optionalNumericFields.every((field) => result[field] === null || isFiniteNumber(result[field])) ||
      (result.energy_evaluations !== null && !Number.isInteger(result.energy_evaluations)) ||
      (result.converged !== null && typeof result.converged !== "boolean") ||
      (result.parameters !== null &&
       (!Array.isArray(result.parameters) || !result.parameters.every(isFiniteNumber)))) {
    throw new Error("The benchmark returned malformed result fields.");
  }
  if (!Object.prototype.hasOwnProperty.call(result, "exact_electronic_energy_hartree") ||
      !Object.prototype.hasOwnProperty.call(result, "reference_electronic_energy_hartree") ||
      !Object.prototype.hasOwnProperty.call(result, "variational_error_hartree") ||
      !Object.prototype.hasOwnProperty.call(result, "energy_evaluations") ||
      !Object.prototype.hasOwnProperty.call(result, "converged") ||
      !Object.prototype.hasOwnProperty.call(result, "parameters")) {
    throw new Error("The benchmark result is missing expected fields.");
  }
  if (Math.abs(result.total_energy_hartree - result.electronic_energy_hartree - result.nuclear_repulsion_hartree) > 1e-10) {
    throw new Error("Total energy does not equal electronic energy plus nuclear repulsion.");
  }
  return result;
}

function addComparisonRow(body, label, value, convention) {
  const row = element("tr");
  row.append(element("td", "", label));
  const numeric = element("td", "text-end font-monospace", value);
  row.append(numeric);
  row.append(element("td", "text-secondary", convention));
  body.append(row);
}

function renderResult(result) {
  document.getElementById("electronic-energy").textContent = formatEnergy(result.electronic_energy_hartree);
  document.getElementById("nuclear-energy").textContent = formatEnergy(result.nuclear_repulsion_hartree);
  document.getElementById("total-energy").textContent = formatEnergy(result.total_energy_hartree);
  document.getElementById("exact-energy").textContent = formatEnergy(result.exact_electronic_energy_hartree);
  document.getElementById("reference-energy").textContent = formatEnergy(result.reference_electronic_energy_hartree);
  document.getElementById("variational-error").textContent = formatError(result.variational_error_hartree);
  document.getElementById("evaluation-count").textContent = formatMaybe(result.energy_evaluations);
  document.getElementById("optimized-parameters").textContent = result.parameters && result.parameters.length
    ? result.parameters.map((value) => value.toPrecision(8)).join(", ")
    : "Unavailable";
  document.getElementById("convergence-detail").textContent =
    result.converged === null ? "Unavailable" : (result.converged ? "Converged" : "Not converged");

  const rows = document.getElementById("comparison-rows");
  rows.replaceChildren();
  addComparisonRow(rows, "Optimized electronic", formatEnergy(result.electronic_energy_hartree), "Variational state · Hartree");
  addComparisonRow(rows, "Exact-sector minimum", formatEnergy(result.exact_electronic_energy_hartree), "Independent small-matrix solve · Hartree");
  addComparisonRow(rows, "Documented electronic reference", formatEnergy(result.reference_electronic_energy_hartree), "Frozen benchmark source · Hartree");
  addComparisonRow(rows, "Nuclear repulsion", formatEnergy(result.nuclear_repulsion_hartree), "Added once to electronic energy");
  addComparisonRow(rows, "Total Born–Oppenheimer", formatEnergy(result.total_energy_hartree), "Electronic + nuclear · Hartree");

  convergenceBadge.className = result.converged === null
    ? "badge bg-secondary-lt"
    : (result.converged ? "badge bg-green-lt" : "badge bg-yellow-lt");
  convergenceBadge.textContent = result.converged === null
    ? "Unavailable"
    : (result.converged ? "Converged" : "Not converged");
  runMessage.hidden = true;
  resultsContent.hidden = false;
}

async function runBenchmark() {
  const id = state.currentId;
  if (!id || !state.benchmarks.has(id)) return;
  const requestId = ++state.requestId;
  resultsContent.hidden = true;
  runMessage.hidden = false;
  runMessage.classList.remove("is-error");
  runMessage.textContent = "Running the QCompile state-vector benchmark…";
  convergenceBadge.className = "badge bg-azure-lt";
  convergenceBadge.textContent = "Running";
  runButton.disabled = true;
  runButton.textContent = "Running…";

  try {
    const response = await fetch(`/api/run?id=${encodeURIComponent(id)}`, { cache: "no-store" });
    const payload = await response.json();
    if (!response.ok) throw new Error(payload && payload.error ? payload.error : `Benchmark request failed (${response.status}).`);
    const result = validateResult(payload, id);
    if (requestId !== state.requestId || id !== state.currentId) return;
    renderResult(result);
  } catch (error) {
    if (requestId !== state.requestId || id !== state.currentId) return;
    resetResults(error instanceof Error ? error.message : "The benchmark failed.", true);
  } finally {
    if (requestId === state.requestId && id === state.currentId) {
      runButton.disabled = false;
      runButton.textContent = "Run again";
    }
  }
}

async function loadCatalog() {
  try {
    const response = await fetch("/api/benchmarks", { cache: "no-store" });
    const payload = await response.json();
    if (!response.ok || !payload || !Array.isArray(payload.benchmarks)) {
      throw new Error(payload && payload.error ? payload.error : "Could not load the benchmark catalog.");
    }
    for (const benchmark of payload.benchmarks) {
      if (!benchmark || typeof benchmark.id !== "string" || state.benchmarks.has(benchmark.id)) {
        throw new Error("The benchmark catalog contains invalid entries.");
      }
      state.benchmarks.set(benchmark.id, benchmark);
      benchmarkList.append(makeBenchmarkCard(benchmark));
    }
    if (state.benchmarks.size !== 2 || !state.benchmarks.has("h2") || !state.benchmarks.has("heh_plus")) {
      throw new Error("Expected the fixed H₂ and HeH⁺ benchmark catalog.");
    }
    catalogStatus.textContent = `${state.benchmarks.size} fixed models`;
  } catch (error) {
    catalogStatus.textContent = "Catalog unavailable";
    benchmarkList.replaceChildren(element("div", "col-12 text-danger", error instanceof Error ? error.message : "Could not load benchmarks."));
    runButton.disabled = true;
  }
}

runButton.addEventListener("click", runBenchmark);
loadCatalog();
