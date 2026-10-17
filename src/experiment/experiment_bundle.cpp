#include "experiment/experiment_bundle.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unistd.h>

namespace {

// JSON Tokenizer & Parser
enum class JsonType { Null, Bool, Number, String, Array, Object };

struct JsonValue {
  JsonType type = JsonType::Null;
  bool boolVal = false;
  double numVal = 0.0;
  std::string strVal;
  std::vector<JsonValue> arrVal;
  std::map<std::string, JsonValue> objVal;

  bool isNull() const { return type == JsonType::Null; }
  bool isBool() const { return type == JsonType::Bool; }
  bool isNumber() const { return type == JsonType::Number; }
  bool isString() const { return type == JsonType::String; }
  bool isArray() const { return type == JsonType::Array; }
  bool isObject() const { return type == JsonType::Object; }

  const JsonValue &operator[](const std::string &key) const {
    static const JsonValue nullVal;
    if (type != JsonType::Object) return nullVal;
    auto it = objVal.find(key);
    return it != objVal.end() ? it->second : nullVal;
  }

  bool contains(const std::string &key) const {
    if (type != JsonType::Object) return false;
    return objVal.find(key) != objVal.end();
  }
};

class JsonParser {
public:
  explicit JsonParser(std::string_view input) : src_(input), pos_(0) {}

  JsonValue parse() {
    skipWhitespace();
    if (pos_ >= src_.size()) {
      throw std::runtime_error("Empty or whitespace-only JSON input");
    }
    JsonValue val = parseValue();
    skipWhitespace();
    if (pos_ < src_.size()) {
      throw std::runtime_error("Unexpected trailing content in JSON input at byte " + std::to_string(pos_));
    }
    return val;
  }

private:
  std::string_view src_;
  std::size_t pos_;

  void skipWhitespace() {
    while (pos_ < src_.size()) {
      const char c = src_[pos_];
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
        ++pos_;
      } else {
        break;
      }
    }
  }

  char peek() const {
    return pos_ < src_.size() ? src_[pos_] : '\0';
  }

  char get() {
    return pos_ < src_.size() ? src_[pos_++] : '\0';
  }

  JsonValue parseValue() {
    skipWhitespace();
    const char c = peek();
    if (c == '{') return parseObject();
    if (c == '[') return parseArray();
    if (c == '"') return parseString();
    if (c == 't' || c == 'f') return parseBool();
    if (c == 'n') return parseNull();
    if (c == '-' || std::isdigit(static_cast<unsigned char>(c))) return parseNumber();
    throw std::runtime_error("Unexpected character '" + std::string(1, c) + "' in JSON at offset " + std::to_string(pos_));
  }

  JsonValue parseObject() {
    get(); // consume '{'
    JsonValue obj;
    obj.type = JsonType::Object;
    skipWhitespace();
    if (peek() == '}') {
      get();
      return obj;
    }
    for (;;) {
      skipWhitespace();
      if (peek() != '"') {
        throw std::runtime_error("Expected string key in JSON object at offset " + std::to_string(pos_));
      }
      JsonValue keyVal = parseString();
      skipWhitespace();
      if (get() != ':') {
        throw std::runtime_error("Expected ':' after key in JSON object at offset " + std::to_string(pos_));
      }
      JsonValue val = parseValue();
      obj.objVal[keyVal.strVal] = std::move(val);
      skipWhitespace();
      const char next = get();
      if (next == '}') break;
      if (next != ',') {
        throw std::runtime_error("Expected ',' or '}' in JSON object at offset " + std::to_string(pos_));
      }
    }
    return obj;
  }

  JsonValue parseArray() {
    get(); // consume '['
    JsonValue arr;
    arr.type = JsonType::Array;
    skipWhitespace();
    if (peek() == ']') {
      get();
      return arr;
    }
    for (;;) {
      JsonValue item = parseValue();
      arr.arrVal.push_back(std::move(item));
      skipWhitespace();
      const char next = get();
      if (next == ']') break;
      if (next != ',') {
        throw std::runtime_error("Expected ',' or ']' in JSON array at offset " + std::to_string(pos_));
      }
    }
    return arr;
  }

  JsonValue parseString() {
    get(); // consume '"'
    std::string str;
    while (pos_ < src_.size()) {
      const char c = get();
      if (c == '"') {
        JsonValue val;
        val.type = JsonType::String;
        val.strVal = std::move(str);
        return val;
      }
      if (c == '\\') {
        const char esc = get();
        switch (esc) {
        case '"': str.push_back('"'); break;
        case '\\': str.push_back('\\'); break;
        case '/': str.push_back('/'); break;
        case 'b': str.push_back('\b'); break;
        case 'f': str.push_back('\f'); break;
        case 'n': str.push_back('\n'); break;
        case 'r': str.push_back('\r'); break;
        case 't': str.push_back('\t'); break;
        case 'u': {
          // 4 hex digits
          std::string hex;
          for (int h = 0; h < 4 && pos_ < src_.size(); ++h) hex.push_back(get());
          unsigned int code = 0;
          std::istringstream hss(hex);
          hss >> std::hex >> code;
          if (code < 128) str.push_back(static_cast<char>(code));
          else str.push_back('?');
          break;
        }
        default: str.push_back(esc); break;
        }
      } else {
        str.push_back(c);
      }
    }
    throw std::runtime_error("Unterminated string in JSON");
  }

  JsonValue parseNumber() {
    std::size_t start = pos_;
    if (peek() == '-') get();
    while (pos_ < src_.size() && std::isdigit(static_cast<unsigned char>(peek()))) get();
    if (peek() == '.') {
      get();
      while (pos_ < src_.size() && std::isdigit(static_cast<unsigned char>(peek()))) get();
    }
    if (peek() == 'e' || peek() == 'E') {
      get();
      if (peek() == '+' || peek() == '-') get();
      while (pos_ < src_.size() && std::isdigit(static_cast<unsigned char>(peek()))) get();
    }
    std::string numStr(src_.substr(start, pos_ - start));
    char *endPtr = nullptr;
    const double val = std::strtod(numStr.c_str(), &endPtr);
    if (!std::isfinite(val)) {
      throw std::runtime_error("JSON numeric literal is not finite: " + numStr);
    }
    JsonValue res;
    res.type = JsonType::Number;
    res.numVal = val;
    return res;
  }

  JsonValue parseBool() {
    if (src_.substr(pos_, 4) == "true") {
      pos_ += 4;
      JsonValue val;
      val.type = JsonType::Bool;
      val.boolVal = true;
      return val;
    }
    if (src_.substr(pos_, 5) == "false") {
      pos_ += 5;
      JsonValue val;
      val.type = JsonType::Bool;
      val.boolVal = false;
      return val;
    }
    throw std::runtime_error("Invalid boolean in JSON at offset " + std::to_string(pos_));
  }

  JsonValue parseNull() {
    if (src_.substr(pos_, 4) == "null") {
      pos_ += 4;
      JsonValue val;
      val.type = JsonType::Null;
      return val;
    }
    throw std::runtime_error("Invalid null literal in JSON at offset " + std::to_string(pos_));
  }
};

void escapeJsonStringTo(std::string_view str, std::ostream &os) {
  os << '"';
  for (char c : str) {
    switch (c) {
    case '"': os << "\\\""; break;
    case '\\': os << "\\\\"; break;
    case '\n': os << "\\n"; break;
    case '\r': os << "\\r"; break;
    case '\t': os << "\\t"; break;
    case '\b': os << "\\b"; break;
    case '\f': os << "\\f"; break;
    default:
      if (static_cast<unsigned char>(c) < 32) {
        os << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c) << std::dec;
      } else {
        os << c;
      }
      break;
    }
  }
  os << '"';
}

std::string escapeJsonString(std::string_view str) {
  std::ostringstream ss;
  escapeJsonStringTo(str, ss);
  return ss.str();
}

} // namespace

std::string serializeExperimentBundle(const ExperimentBundle &b) {
  std::ostringstream os;
  os << std::setprecision(17);
  os << "{\n";
  os << "  \"schema\": " << escapeJsonString(b.schema) << ",\n";
  os << "  \"metadata\": {\n";
  os << "    \"id\": " << escapeJsonString(b.id) << ",\n";
  os << "    \"title\": " << escapeJsonString(b.title) << ",\n";
  os << "    \"timestamp\": " << escapeJsonString(b.timestamp) << ",\n";
  os << "    \"description\": " << escapeJsonString(b.description) << "\n";
  os << "  },\n";

  // Source & compiler config
  os << "  \"source\": {\n";
  os << "    \"dsl\": " << escapeJsonString(b.dslSource) << ",\n";
  os << "    \"optimize\": " << (b.optimize ? "true" : "false") << ",\n";
  os << "    \"diagnostics\": [";
  for (std::size_t i = 0; i < b.diagnostics.size(); ++i) {
    if (i != 0) os << ", ";
    escapeJsonStringTo(b.diagnostics[i], os);
  }
  os << "]\n  },\n";

  // Circuit model
  os << "  \"circuit\": {\n";
  os << "    \"qubit_count\": " << b.qubitCount << ",\n";
  auto writeInstrArray = [&](const std::vector<Instruction> &instrs) {
    os << "[\n";
    for (std::size_t i = 0; i < instrs.size(); ++i) {
      const auto &in = instrs[i];
      os << "      {\"type\": " << (in.type == InstructionType::Measure ? "\"Measure\"" : "\"ApplyGate\"")
         << ", \"gate\": " << escapeJsonString(in.gateName)
         << ", \"qubits\": [";
      for (std::size_t q = 0; q < in.qubitIndices.size(); ++q) {
        if (q != 0) os << ", ";
        os << in.qubitIndices[q];
      }
      os << "], \"parameter\": " << in.parameterValue << "}";
      if (i + 1 < instrs.size()) os << ",";
      os << "\n";
    }
    os << "    ]";
  };
  os << "    \"raw_instructions\": ";
  writeInstrArray(b.rawInstructions);
  os << ",\n";
  os << "    \"instructions\": ";
  writeInstrArray(b.instructions);
  os << "\n  },\n";

  // Simulation state
  os << "  \"simulation\": {\n";
  os << "    \"current_step\": " << b.currentStep << ",\n";
  os << "    \"total_steps\": " << b.totalSteps << ",\n";
  os << "    \"is_complete\": " << (b.isComplete ? "true" : "false") << ",\n";
  os << "    \"state_norm\": " << b.stateNorm << ",\n";
  os << "    \"statevector\": [";
  for (std::size_t i = 0; i < b.statevector.size(); ++i) {
    if (i != 0) os << ", ";
    os << "{\"re\": " << b.statevector[i].real() << ", \"im\": " << b.statevector[i].imag() << "}";
  }
  os << "],\n";
  os << "    \"basis_probabilities\": [";
  for (std::size_t i = 0; i < b.basisProbabilities.size(); ++i) {
    if (i != 0) os << ", ";
    os << b.basisProbabilities[i];
  }
  os << "],\n";
  os << "    \"recorded_measurements\": [";
  for (std::size_t i = 0; i < b.recordedMeasurements.size(); ++i) {
    if (i != 0) os << ", ";
    os << "{\"step\": " << b.recordedMeasurements[i].first << ", \"outcome\": "
       << escapeJsonString(b.recordedMeasurements[i].second) << "}";
  }
  os << "]\n  },\n";

  // Chemistry
  os << "  \"chemistry\": {\n";
  os << "    \"has_chemistry\": " << (b.hasChemistry ? "true" : "false");
  if (b.hasChemistry) {
    os << ",\n    \"molecule_name\": " << escapeJsonString(b.moleculeName)
       << ",\n    \"formula\": " << escapeJsonString(b.formula)
       << ",\n    \"basis\": " << escapeJsonString(b.basis)
       << ",\n    \"mapping\": " << escapeJsonString(b.mapping)
       << ",\n    \"distance\": " << b.distance
       << ",\n    \"distance_units\": " << escapeJsonString(b.distanceUnits)
       << ",\n    \"electronic_energy\": " << b.electronicEnergy
       << ",\n    \"nuclear_repulsion\": " << b.nuclearRepulsion
       << ",\n    \"total_energy\": " << b.totalEnergy
       << ",\n    \"exact_electronic_energy\": " << b.exactElectronicEnergy
       << ",\n    \"exact_total_energy\": " << b.exactTotalEnergy
       << ",\n    \"variational_error\": " << b.variationalError
       << ",\n    \"evaluations\": " << b.energyEvaluations
       << ",\n    \"converged\": " << (b.converged ? "true" : "false");
  }
  os << "\n  },\n";

  // Convergence trajectory
  os << "  \"convergence\": {\n";
  os << "    \"has_convergence\": " << (b.hasConvergence ? "true" : "false");
  if (b.hasConvergence) {
    os << ",\n    \"trajectory\": [\n";
    for (std::size_t i = 0; i < b.convergenceTrajectory.size(); ++i) {
      const auto &p = b.convergenceTrajectory[i];
      os << "      {\"eval\": " << p.evaluation << ", \"p0\": " << p.parameter0 << ", \"p1\": " << p.parameter1
         << ", \"electronic\": " << p.electronicEnergy << ", \"total\": " << p.totalEnergy
         << ", \"is_best\": " << (p.isBest ? "true" : "false") << "}";
      if (i + 1 < b.convergenceTrajectory.size()) os << ",";
      os << "\n";
    }
    os << "    ]";
  }
  os << "\n  }\n";
  os << "}\n";
  return os.str();
}

ExperimentBundle deserializeExperimentBundle(const std::string &jsonText) {
  JsonParser parser(jsonText);
  JsonValue root = parser.parse();

  if (!root.isObject()) {
    throw std::runtime_error("Invalid .qcomp experiment bundle: root must be a JSON object");
  }

  // Schema validation
  if (!root.contains("schema") || !root["schema"].isString()) {
    throw std::runtime_error("Malformed .qcomp bundle: missing or invalid 'schema' identifier");
  }
  const std::string schema = root["schema"].strVal;
  if (schema != "QCOMPILE-EXPERIMENT-1") {
    throw std::runtime_error("Unsupported .qcomp bundle schema version: '" + schema + "' (expected QCOMPILE-EXPERIMENT-1)");
  }

  ExperimentBundle b;
  b.schema = schema;

  // Metadata
  const auto &meta = root["metadata"];
  if (meta.isObject()) {
    if (meta.contains("id") && meta["id"].isString()) b.id = meta["id"].strVal;
    if (meta.contains("title") && meta["title"].isString()) b.title = meta["title"].strVal;
    if (meta.contains("timestamp") && meta["timestamp"].isString()) b.timestamp = meta["timestamp"].strVal;
    if (meta.contains("description") && meta["description"].isString()) b.description = meta["description"].strVal;
  }

  // Source
  const auto &source = root["source"];
  if (!source.isObject() || !source.contains("dsl") || !source["dsl"].isString()) {
    throw std::runtime_error("Malformed .qcomp bundle: missing or invalid 'source.dsl' field");
  }
  b.dslSource = source["dsl"].strVal;
  if (source.contains("optimize") && source["optimize"].isBool()) {
    b.optimize = source["optimize"].boolVal;
  }
  if (source.contains("diagnostics") && source["diagnostics"].isArray()) {
    for (const auto &item : source["diagnostics"].arrVal) {
      if (item.isString()) b.diagnostics.push_back(item.strVal);
    }
  }

  // Circuit
  const auto &circuit = root["circuit"];
  if (!circuit.isObject() || !circuit.contains("qubit_count") || !circuit["qubit_count"].isNumber()) {
    throw std::runtime_error("Malformed .qcomp bundle: missing or invalid 'circuit.qubit_count' field");
  }
  b.qubitCount = static_cast<int>(circuit["qubit_count"].numVal);
  if (b.qubitCount < 0 || b.qubitCount > 24) {
    throw std::runtime_error("Invalid qubit count in .qcomp bundle: " + std::to_string(b.qubitCount));
  }

  auto parseInstrVector = [&](const JsonValue &arrNode, std::vector<Instruction> &out) {
    if (!arrNode.isArray()) return;
    for (const auto &item : arrNode.arrVal) {
      if (!item.isObject()) continue;
      Instruction in;
      std::string t = item["type"].isString() ? item["type"].strVal : "ApplyGate";
      in.type = (t == "Measure") ? InstructionType::Measure : InstructionType::ApplyGate;
      in.gateName = item["gate"].isString() ? item["gate"].strVal : "";
      if (item.contains("parameter") && item["parameter"].isNumber()) {
        in.parameterValue = item["parameter"].numVal;
      }
      if (item.contains("qubits") && item["qubits"].isArray()) {
        for (const auto &q : item["qubits"].arrVal) {
          if (q.isNumber()) in.qubitIndices.push_back(static_cast<int>(q.numVal));
        }
      }
      out.push_back(std::move(in));
    }
  };

  parseInstrVector(circuit["raw_instructions"], b.rawInstructions);
  parseInstrVector(circuit["instructions"], b.instructions);

  // Simulation
  const auto &sim = root["simulation"];
  if (sim.isObject()) {
    if (sim.contains("current_step") && sim["current_step"].isNumber()) {
      b.currentStep = static_cast<std::size_t>(sim["current_step"].numVal);
    }
    if (sim.contains("total_steps") && sim["total_steps"].isNumber()) {
      b.totalSteps = static_cast<std::size_t>(sim["total_steps"].numVal);
    }
    if (sim.contains("is_complete") && sim["is_complete"].isBool()) {
      b.isComplete = sim["is_complete"].boolVal;
    }
    if (sim.contains("state_norm") && sim["state_norm"].isNumber()) {
      b.stateNorm = sim["state_norm"].numVal;
    }
    if (sim.contains("statevector") && sim["statevector"].isArray()) {
      for (const auto &item : sim["statevector"].arrVal) {
        if (item.isObject()) {
          double re = item["re"].isNumber() ? item["re"].numVal : 0.0;
          double im = item["im"].isNumber() ? item["im"].numVal : 0.0;
          b.statevector.push_back(std::complex<double>(re, im));
        }
      }
    }
    if (sim.contains("basis_probabilities") && sim["basis_probabilities"].isArray()) {
      for (const auto &item : sim["basis_probabilities"].arrVal) {
        if (item.isNumber()) b.basisProbabilities.push_back(item.numVal);
      }
    }
    if (sim.contains("recorded_measurements") && sim["recorded_measurements"].isArray()) {
      for (const auto &item : sim["recorded_measurements"].arrVal) {
        if (item.isObject() && item.contains("step") && item.contains("outcome")) {
          std::size_t step = static_cast<std::size_t>(item["step"].numVal);
          std::string outcome = item["outcome"].isString() ? item["outcome"].strVal : "";
          b.recordedMeasurements.emplace_back(step, outcome);
        }
      }
    }
  }

  // Chemistry
  const auto &chem = root["chemistry"];
  if (chem.isObject() && chem.contains("has_chemistry") && chem["has_chemistry"].isBool()) {
    b.hasChemistry = chem["has_chemistry"].boolVal;
    if (b.hasChemistry) {
      if (chem.contains("molecule_name") && chem["molecule_name"].isString()) b.moleculeName = chem["molecule_name"].strVal;
      if (chem.contains("formula") && chem["formula"].isString()) b.formula = chem["formula"].strVal;
      if (chem.contains("basis") && chem["basis"].isString()) b.basis = chem["basis"].strVal;
      if (chem.contains("mapping") && chem["mapping"].isString()) b.mapping = chem["mapping"].strVal;
      if (chem.contains("distance") && chem["distance"].isNumber()) b.distance = chem["distance"].numVal;
      if (chem.contains("distance_units") && chem["distance_units"].isString()) b.distanceUnits = chem["distance_units"].strVal;
      if (chem.contains("electronic_energy") && chem["electronic_energy"].isNumber()) b.electronicEnergy = chem["electronic_energy"].numVal;
      if (chem.contains("nuclear_repulsion") && chem["nuclear_repulsion"].isNumber()) b.nuclearRepulsion = chem["nuclear_repulsion"].numVal;
      if (chem.contains("total_energy") && chem["total_energy"].isNumber()) b.totalEnergy = chem["total_energy"].numVal;
      if (chem.contains("exact_electronic_energy") && chem["exact_electronic_energy"].isNumber()) b.exactElectronicEnergy = chem["exact_electronic_energy"].numVal;
      if (chem.contains("exact_total_energy") && chem["exact_total_energy"].isNumber()) b.exactTotalEnergy = chem["exact_total_energy"].numVal;
      if (chem.contains("variational_error") && chem["variational_error"].isNumber()) b.variationalError = chem["variational_error"].numVal;
      if (chem.contains("evaluations") && chem["evaluations"].isNumber()) b.energyEvaluations = static_cast<int>(chem["evaluations"].numVal);
      if (chem.contains("converged") && chem["converged"].isBool()) b.converged = chem["converged"].boolVal;
    }
  }

  // Convergence
  const auto &conv = root["convergence"];
  if (conv.isObject() && conv.contains("has_convergence") && conv["has_convergence"].isBool()) {
    b.hasConvergence = conv["has_convergence"].boolVal;
    if (b.hasConvergence && conv.contains("trajectory") && conv["trajectory"].isArray()) {
      for (const auto &item : conv["trajectory"].arrVal) {
        if (item.isObject()) {
          VqeConvergencePoint pt;
          if (item.contains("eval") && item["eval"].isNumber()) pt.evaluation = static_cast<int>(item["eval"].numVal);
          if (item.contains("p0") && item["p0"].isNumber()) pt.parameter0 = item["p0"].numVal;
          if (item.contains("p1") && item["p1"].isNumber()) pt.parameter1 = item["p1"].numVal;
          if (item.contains("electronic") && item["electronic"].isNumber()) pt.electronicEnergy = item["electronic"].numVal;
          if (item.contains("total") && item["total"].isNumber()) pt.totalEnergy = item["total"].numVal;
          if (item.contains("is_best") && item["is_best"].isBool()) pt.isBest = item["is_best"].boolVal;
          b.convergenceTrajectory.push_back(pt);
        }
      }
    }
  }

  return b;
}

void saveExperimentBundle(const std::string &filePath, const ExperimentBundle &bundle) {
  static std::atomic<unsigned long> tempIdCounter{0};
  const std::string tempPath = filePath + ".tmp." + std::to_string(getpid()) + "." + std::to_string(tempIdCounter.fetch_add(1));

  const std::string content = serializeExperimentBundle(bundle);
  {
    std::ofstream out(tempPath, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
      throw std::runtime_error("Could not open temporary file for writing: " + tempPath);
    }
    out << content;
    if (!out) {
      out.close();
      std::remove(tempPath.c_str());
      throw std::runtime_error("Failed to write experiment bundle content to temporary file: " + tempPath);
    }
  }

  // Atomic file replace
  std::error_code ec;
  std::filesystem::rename(tempPath, filePath, ec);
  if (ec) {
    // Fallback if destination exists and cannot be overwritten directly by rename
    std::filesystem::remove(filePath, ec);
    std::filesystem::rename(tempPath, filePath, ec);
    if (ec) {
      std::remove(tempPath.c_str());
      throw std::runtime_error("Failed to atomically replace destination bundle: " + filePath + " (" + ec.message() + ")");
    }
  }
}

ExperimentBundle loadExperimentBundle(const std::string &filePath) {
  std::ifstream in(filePath, std::ios::in);
  if (!in.is_open()) {
    throw std::runtime_error("Could not open .qcomp bundle file: " + filePath);
  }
  std::ostringstream ss;
  ss << in.rdbuf();
  if (!in && !in.eof()) {
    throw std::runtime_error("Error reading .qcomp bundle file: " + filePath);
  }
  return deserializeExperimentBundle(ss.str());
}

std::string exportConvergenceCsv(const ExperimentBundle &bundle) {
  std::ostringstream ss;
  ss << "evaluation,parameter0,parameter1,electronic_energy_hartree,nuclear_repulsion_hartree,total_energy_hartree,is_best_so_far\n";
  ss << std::setprecision(17);

  const double nucRep = bundle.hasChemistry ? bundle.nuclearRepulsion : 0.0;
  for (const auto &p : bundle.convergenceTrajectory) {
    ss << p.evaluation << ","
       << p.parameter0 << ","
       << p.parameter1 << ","
       << p.electronicEnergy << ","
       << nucRep << ","
       << p.totalEnergy << ","
       << (p.isBest ? "true" : "false") << "\n";
  }
  return ss.str();
}

void saveConvergenceCsv(const std::string &filePath, const ExperimentBundle &bundle) {
  static std::atomic<unsigned long> tempIdCounter{0};
  const std::string tempPath = filePath + ".tmp." + std::to_string(getpid()) + "." + std::to_string(tempIdCounter.fetch_add(1));

  const std::string content = exportConvergenceCsv(bundle);
  {
    std::ofstream out(tempPath, std::ios::out | std::ios::trunc);
    if (!out.is_open()) throw std::runtime_error("Could not create temporary CSV: " + tempPath);
    out << content;
    if (!out) {
      out.close();
      std::remove(tempPath.c_str());
      throw std::runtime_error("Failed writing CSV content to: " + tempPath);
    }
  }

  std::error_code ec;
  std::filesystem::rename(tempPath, filePath, ec);
  if (ec) {
    std::filesystem::remove(filePath, ec);
    std::filesystem::rename(tempPath, filePath, ec);
    if (ec) {
      std::remove(tempPath.c_str());
      throw std::runtime_error("Failed replacing CSV file: " + filePath + " (" + ec.message() + ")");
    }
  }
}

std::string exportCircuitSvg(const std::vector<Instruction> &instructions, int qubitCount) {
  constexpr int kWireSpacing = 52;
  constexpr int kLeftMargin = 64;
  constexpr int kTopMargin = 42;
  constexpr int kColumnSpacing = 68;

  const int fullWidth = std::max(520, kLeftMargin + static_cast<int>(instructions.size() + 1) * kColumnSpacing + 48);
  const int fullHeight = std::max(180, kTopMargin + std::max(1, qubitCount) * kWireSpacing + 28);

  std::ostringstream os;
  os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << fullWidth << " " << fullHeight
     << "\" width=\"" << fullWidth << "\" height=\"" << fullHeight << "\">\n";

  // Technical dark canvas background
  os << "  <rect width=\"100%\" height=\"100%\" fill=\"#0b0f17\"/>\n";

  if (qubitCount <= 0) {
    os << "  <text x=\"" << (fullWidth / 2) << "\" y=\"" << (fullHeight / 2)
       << "\" fill=\"#64748b\" font-family=\"monospace\" font-size=\"12\" text-anchor=\"middle\">"
       << "No quantum register declared</text>\n";
    os << "</svg>\n";
    return os.str();
  }

  // Wires & labels
  for (int q = 0; q < qubitCount; ++q) {
    const int y = kTopMargin + q * kWireSpacing;
    os << "  <!-- Qubit q" << q << " -->\n";
    os << "  <rect x=\"12\" y=\"" << (y - 10) << "\" width=\"36\" height=\"20\" rx=\"2\" fill=\"#0f172a\" stroke=\"#334155\" stroke-width=\"1\"/>\n";
    os << "  <text x=\"30\" y=\"" << (y + 4) << "\" fill=\"#94a3b8\" font-family=\"monospace\" font-size=\"10\" font-weight=\"bold\" text-anchor=\"middle\">q" << q << "</text>\n";
    os << "  <line x1=\"" << (kLeftMargin - 12) << "\" y1=\"" << y << "\" x2=\"" << (fullWidth - 20) << "\" y2=\"" << y << "\" stroke=\"#334155\" stroke-width=\"2\"/>\n";
  }

  // Gate instructions
  for (std::size_t i = 0; i < instructions.size(); ++i) {
    const auto &in = instructions[i];
    const int x = kLeftMargin + static_cast<int>(i) * kColumnSpacing + 32;

    if (in.type == InstructionType::Measure) {
      for (int q : in.qubitIndices) {
        if (q < 0 || q >= qubitCount) continue;
        const int y = kTopMargin + q * kWireSpacing;
        os << "  <rect x=\"" << (x - 20) << "\" y=\"" << (y - 14) << "\" width=\"40\" height=\"28\" rx=\"4\" fill=\"#1e1b4b\" stroke=\"#818cf8\" stroke-width=\"1.5\"/>\n";
        os << "  <text x=\"" << x << "\" y=\"" << (y + 4) << "\" fill=\"#c7d2fe\" font-family=\"monospace\" font-size=\"9\" font-weight=\"bold\" text-anchor=\"middle\">M</text>\n";
      }
    } else if (in.gateName == "CNOT" && in.qubitIndices.size() >= 2) {
      const int cQ = in.qubitIndices[0];
      const int tQ = in.qubitIndices[1];
      if (cQ >= 0 && cQ < qubitCount && tQ >= 0 && tQ < qubitCount) {
        const int cY = kTopMargin + cQ * kWireSpacing;
        const int tY = kTopMargin + tQ * kWireSpacing;
        // Bus
        os << "  <line x1=\"" << x << "\" y1=\"" << cY << "\" x2=\"" << x << "\" y2=\"" << tY << "\" stroke=\"#38bdf8\" stroke-width=\"2\"/>\n";
        // Control dot
        os << "  <circle cx=\"" << x << "\" cy=\"" << cY << "\" r=\"5\" fill=\"#38bdf8\"/>\n";
        // Target circle with cross
        os << "  <circle cx=\"" << x << "\" cy=\"" << tY << "\" r=\"10\" fill=\"#0f172a\" stroke=\"#38bdf8\" stroke-width=\"1.5\"/>\n";
        os << "  <line x1=\"" << (x - 10) << "\" y1=\"" << tY << "\" x2=\"" << (x + 10) << "\" y2=\"" << tY << "\" stroke=\"#38bdf8\" stroke-width=\"1.5\"/>\n";
        os << "  <line x1=\"" << x << "\" y1=\"" << (tY - 10) << "\" x2=\"" << x << "\" y2=\"" << (tY + 10) << "\" stroke=\"#38bdf8\" stroke-width=\"1.5\"/>\n";
      }
    } else if (!in.qubitIndices.empty()) {
      const int q = in.qubitIndices[0];
      if (q >= 0 && q < qubitCount) {
        const int y = kTopMargin + q * kWireSpacing;
        std::string label = in.gateName;
        if (in.gateName == "RZ") {
          std::ostringstream ss;
          ss << "RZ(" << std::fixed << std::setprecision(2) << in.parameterValue << ")";
          label = ss.str();
        }
        const int boxW = std::max(40, static_cast<int>(label.size() * 8 + 12));
        os << "  <rect x=\"" << (x - boxW / 2) << "\" y=\"" << (y - 14) << "\" width=\"" << boxW << "\" height=\"28\" rx=\"4\" fill=\"#0f172a\" stroke=\"#38bdf8\" stroke-width=\"1.5\"/>\n";
        os << "  <text x=\"" << x << "\" y=\"" << (y + 4) << "\" fill=\"#f1f5f9\" font-family=\"monospace\" font-size=\"9\" font-weight=\"bold\" text-anchor=\"middle\">" << label << "</text>\n";
      }
    }
  }

  os << "</svg>\n";
  return os.str();
}

void saveCircuitSvg(const std::string &filePath, const std::vector<Instruction> &instructions, int qubitCount) {
  static std::atomic<unsigned long> tempIdCounter{0};
  const std::string tempPath = filePath + ".tmp." + std::to_string(getpid()) + "." + std::to_string(tempIdCounter.fetch_add(1));

  const std::string svg = exportCircuitSvg(instructions, qubitCount);
  {
    std::ofstream out(tempPath, std::ios::out | std::ios::trunc);
    if (!out.is_open()) throw std::runtime_error("Could not create temporary SVG file: " + tempPath);
    out << svg;
    if (!out) {
      out.close();
      std::remove(tempPath.c_str());
      throw std::runtime_error("Failed writing SVG content to: " + tempPath);
    }
  }

  std::error_code ec;
  std::filesystem::rename(tempPath, filePath, ec);
  if (ec) {
    std::filesystem::remove(filePath, ec);
    std::filesystem::rename(tempPath, filePath, ec);
    if (ec) {
      std::remove(tempPath.c_str());
      throw std::runtime_error("Failed replacing SVG file: " + filePath + " (" + ec.message() + ")");
    }
  }
}

ExperimentComparisonResult compareExperiments(const ExperimentBundle &a, const ExperimentBundle &b) {
  ExperimentComparisonResult res;
  res.titleA = a.title.empty() ? (a.id.empty() ? "Experiment A" : a.id) : a.title;
  res.titleB = b.title.empty() ? (b.id.empty() ? "Experiment B" : b.id) : b.title;
  res.formulaA = a.hasChemistry ? a.formula : "Pure Circuit";
  res.formulaB = b.hasChemistry ? b.formula : "Pure Circuit";

  res.moleculesMatch = (a.hasChemistry == b.hasChemistry) && (!a.hasChemistry || a.formula == b.formula);
  res.geometriesMatch = (a.hasChemistry == b.hasChemistry) && (!a.hasChemistry || (std::abs(a.distance - b.distance) < 1e-6 && a.distanceUnits == b.distanceUnits));

  if (!a.hasChemistry && !b.hasChemistry) {
    res.isDirectlyComparable = true;
    res.summary = "Both experiments are pure quantum circuit simulations.";
  } else if (a.hasChemistry && b.hasChemistry) {
    if (a.formula != b.formula) {
      res.isDirectlyComparable = false;
      res.summary = "Mismatched molecules (" + a.formula + " vs " + b.formula + "): energies reflect distinct chemical systems.";
    } else if (std::abs(a.distance - b.distance) > 1e-4) {
      res.isDirectlyComparable = true;
      res.summary = "Same molecule (" + a.formula + ") at distinct bond lengths: comparing potential energy surface points.";
    } else {
      res.isDirectlyComparable = true;
      res.summary = "Same molecule (" + a.formula + ") at identical geometry: directly comparable VQE runs.";
    }
  } else {
    res.isDirectlyComparable = false;
    res.summary = "Incompatible run types: one is a pure circuit simulation and the other is a molecular chemistry experiment.";
  }

  auto addComparison = [&](std::string_view prop, const std::string &valA, const std::string &valB,
                           const std::string &delta, bool compatible, std::string_view note) {
    res.items.push_back({std::string(prop), valA, valB, delta, compatible, std::string(note)});
  };

  // 1. Qubit & Circuit Metrics
  addComparison("Qubit Count", std::to_string(a.qubitCount), std::to_string(b.qubitCount),
                std::to_string(b.qubitCount - a.qubitCount), a.qubitCount == b.qubitCount,
                a.qubitCount == b.qubitCount ? "Identical register size" : "Different register dimensions");

  addComparison("Instruction Count (Optimized)", std::to_string(a.instructions.size()), std::to_string(b.instructions.size()),
                std::to_string(static_cast<long>(b.instructions.size()) - static_cast<long>(a.instructions.size())),
                true, "Circuit length comparison");

  addComparison("State Norm", std::to_string(a.stateNorm), std::to_string(b.stateNorm),
                std::to_string(b.stateNorm - a.stateNorm),
                std::abs(a.stateNorm - 1.0) < 1e-10 && std::abs(b.stateNorm - 1.0) < 1e-10,
                "Both states normalized");

  // 2. Chemistry Metrics if present
  if (a.hasChemistry || b.hasChemistry) {
    addComparison("Molecule", a.hasChemistry ? a.formula : "N/A", b.hasChemistry ? b.formula : "N/A",
                  "—", res.moleculesMatch, res.moleculesMatch ? "Same chemical species" : "Mismatched species");

    addComparison("Basis Set", a.hasChemistry ? a.basis : "N/A", b.hasChemistry ? b.basis : "N/A",
                  "—", a.hasChemistry && b.hasChemistry && a.basis == b.basis, "Orbital basis comparison");

    std::ostringstream dA, dB, dDelta;
    dA << std::setprecision(6) << a.distance << " " << a.distanceUnits;
    dB << std::setprecision(6) << b.distance << " " << b.distanceUnits;
    dDelta << (b.distance - a.distance) << " " << a.distanceUnits;
    addComparison("Bond Distance", a.hasChemistry ? dA.str() : "N/A", b.hasChemistry ? dB.str() : "N/A",
                  (a.hasChemistry && b.hasChemistry && a.distanceUnits == b.distanceUnits) ? dDelta.str() : "—",
                  res.geometriesMatch, res.geometriesMatch ? "Identical coordinates" : "Geometry delta");

    std::ostringstream eTotA, eTotB, dEtot;
    eTotA << std::setprecision(10) << a.totalEnergy << " Ha";
    eTotB << std::setprecision(10) << b.totalEnergy << " Ha";
    dEtot << std::setprecision(10) << (b.totalEnergy - a.totalEnergy) << " Ha";
    addComparison("Total Energy (Born–Oppenheimer)", a.hasChemistry ? eTotA.str() : "N/A", b.hasChemistry ? eTotB.str() : "N/A",
                  (a.hasChemistry && b.hasChemistry) ? dEtot.str() : "—",
                  res.isDirectlyComparable, "Total energy comparison");

    std::ostringstream varA, varB;
    varA << std::setprecision(6) << a.variationalError << " Ha";
    varB << std::setprecision(6) << b.variationalError << " Ha";
    addComparison("Variational Error (ΔE)", a.hasChemistry ? varA.str() : "N/A", b.hasChemistry ? varB.str() : "N/A",
                  "—", true, "Error against physical-sector exact ground state");

    addComparison("Optimizer Evaluations", a.hasChemistry ? std::to_string(a.energyEvaluations) : "N/A",
                  b.hasChemistry ? std::to_string(b.energyEvaluations) : "N/A",
                  (a.hasChemistry && b.hasChemistry) ? std::to_string(b.energyEvaluations - a.energyEvaluations) : "—",
                  true, "Evaluation count to reach convergence");
  }

  return res;
}
