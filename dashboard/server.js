"use strict";

const childProcess = require("node:child_process");
const fs = require("node:fs");
const http = require("node:http");
const path = require("node:path");

const projectRoot = path.resolve(__dirname, "..");
const executable = path.join(projectRoot, "qcompile_benchmarks");
const allowedIds = new Set(["h2", "heh_plus"]);
const staticFiles = new Map([
  ["/", [path.join(__dirname, "index.html"), "text/html; charset=utf-8"]],
  ["/app.css", [path.join(__dirname, "app.css"), "text/css; charset=utf-8"]],
  ["/app.js", [path.join(__dirname, "app.js"), "text/javascript; charset=utf-8"]],
  ["/vendor/tabler.min.css", [path.join(__dirname, "vendor", "tabler.min.css"), "text/css; charset=utf-8"]]
]);

function send(response, statusCode, contentType, body) {
  response.writeHead(statusCode, {
    "Content-Type": contentType,
    "Cache-Control": "no-store",
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
    "Content-Security-Policy": "default-src 'self'; script-src 'self'; style-src 'self'; connect-src 'self'; img-src 'self' data:; object-src 'none'; base-uri 'none'; frame-ancestors 'none'"
  });
  response.end(body);
}

function sendJson(response, statusCode, value) {
  send(response, statusCode, "application/json; charset=utf-8", JSON.stringify(value));
}

function runExecutable(args) {
  return new Promise((resolve, reject) => {
    childProcess.execFile(executable, args, {
      cwd: projectRoot,
      encoding: "utf8",
      timeout: 30000,
      maxBuffer: 256 * 1024,
      windowsHide: true,
      shell: false
    }, (error, stdout, stderr) => {
      if (error) {
        const detail = (stderr || error.message || "Benchmark process failed.").trim();
        reject(new Error(detail));
        return;
      }
      resolve(stdout);
    });
  });
}

async function handleRequest(request, response) {
  if (request.method !== "GET") {
    sendJson(response, 405, { error: "Only GET requests are supported." });
    return;
  }

  const url = new URL(request.url, "http://127.0.0.1");
  if (url.pathname === "/api/benchmarks") {
    try {
      const catalog = JSON.parse(await runExecutable(["--json-list"]));
      if (!catalog || !Array.isArray(catalog.benchmarks) ||
          catalog.benchmarks.length !== allowedIds.size ||
          catalog.benchmarks.some((benchmark) => !benchmark || !allowedIds.has(benchmark.id))) {
        throw new Error("The benchmark executable returned an invalid catalog.");
      }
      sendJson(response, 200, catalog);
    } catch (error) {
      sendJson(response, 502, { error: error instanceof Error ? error.message : "Could not read benchmark catalog." });
    }
    return;
  }

  if (url.pathname === "/api/run") {
    const id = url.searchParams.get("id");
    if (!allowedIds.has(id)) {
      sendJson(response, 400, { error: "Unknown benchmark selection." });
      return;
    }
    try {
      const result = JSON.parse(await runExecutable(["--json", id]));
      if (!result || !result.benchmark || result.benchmark.id !== id) {
        throw new Error("The benchmark executable returned malformed output.");
      }
      sendJson(response, 200, result);
    } catch (error) {
      sendJson(response, 502, { error: error instanceof Error ? error.message : "Benchmark execution failed." });
    }
    return;
  }

  const asset = staticFiles.get(url.pathname);
  if (!asset) {
    send(response, 404, "text/plain; charset=utf-8", "Not found");
    return;
  }
  fs.readFile(asset[0], (error, content) => {
    if (error) {
      send(response, 500, "text/plain; charset=utf-8", "Dashboard asset unavailable");
      return;
    }
    send(response, 200, asset[1], content);
  });
}

const configuredPort = Number(process.env.PORT || 4173);
if (!Number.isInteger(configuredPort) || configuredPort < 1024 || configuredPort > 65535) {
  console.error("PORT must be an integer from 1024 through 65535.");
  process.exit(2);
}

const server = http.createServer((request, response) => {
  handleRequest(request, response).catch((error) => {
    if (!response.headersSent) {
      sendJson(response, 500, { error: error instanceof Error ? error.message : "Unexpected server error." });
    } else {
      response.destroy();
    }
  });
});

server.on("error", (error) => {
  console.error(`Dashboard server failed: ${error.message}`);
  process.exitCode = 1;
});

server.listen(configuredPort, "127.0.0.1", () => {
  console.log(`QCompile dashboard listening at http://127.0.0.1:${configuredPort}`);
});
