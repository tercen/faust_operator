//! Operator properties. The three the R operator has (`selectionQuantile`,
//! `depthScoreThreshold`, `seed`) with the same names and defaults, plus `sample_factor` (the
//! R operator hard-codes `filename`) and `threads`.
use anyhow::{Result, bail};
use tercen_rs::PropertyReader;
use tercen_rs::context::ContextBase;

#[derive(Debug, Clone)]
pub struct Settings {
    /// `selectionQuantile`: a marker is kept when this quantile of its depth score across
    /// samples is at least `depthScoreThreshold`.
    pub selection_quantile: f64,
    /// `depthScoreThreshold`.
    pub depth_score_threshold: f64,
    /// `seed`: `seedValue` of FAUST — the forest seed is derived from it the way R does.
    pub seed: u32,
    /// `sample_factor`: the column factor that names the sample. Empty means `filename` if
    /// projected, else one sample.
    pub sample_factor: String,
    /// `threads`: samples are processed in parallel; 0 uses every core the task has.
    pub threads: usize,
    /// `histograms`: draw the per-marker threshold histograms (the R operator's third table).
    pub histograms: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            selection_quantile: 0.5,
            depth_score_threshold: 0.01,
            seed: 123,
            sample_factor: String::new(),
            threads: 0,
            histograms: true,
        }
    }
}

pub fn read(ctx: &ContextBase) -> Result<Settings> {
    let pr = PropertyReader::from_operator_settings(ctx.operator_settings());
    let d = Settings::default();
    let num = |name: &str, dflt: f64| -> Result<f64> {
        let raw = pr.get_string(name, "").trim().to_string();
        if raw.is_empty() || raw.eq_ignore_ascii_case("null") {
            return Ok(dflt);
        }
        raw.parse::<f64>()
            .map_err(|_| anyhow::anyhow!("property '{name}' should be a number, got '{raw}'"))
    };
    let text = |name: &str| pr.get_string(name, "").trim().to_string();
    let selection_quantile = num("selectionQuantile", d.selection_quantile)?;
    if !(0.0..=1.0).contains(&selection_quantile) {
        bail!("property 'selectionQuantile' must be between 0 and 1, got {selection_quantile}");
    }
    let depth_score_threshold = num("depthScoreThreshold", d.depth_score_threshold)?;
    if !(0.0..=1.0).contains(&depth_score_threshold) {
        bail!(
            "property 'depthScoreThreshold' must be between 0 and 1, got {depth_score_threshold}"
        );
    }
    let seed = num("seed", d.seed as f64)?;
    if !(1.0..=4_294_967_295.0).contains(&seed) {
        bail!("property 'seed' must be a positive integer below 2^32, got {seed}");
    }
    let histograms = match text("histograms").to_lowercase().as_str() {
        "" => d.histograms,
        "true" | "yes" | "1" => true,
        "false" | "no" | "0" => false,
        other => bail!("property 'histograms' must be true or false, got '{other}'"),
    };
    let s = Settings {
        selection_quantile,
        depth_score_threshold,
        seed: seed as u32,
        sample_factor: text("sample_factor"),
        threads: num("threads", 0.0)?.max(0.0) as usize,
        histograms,
    };
    tracing::info!(?s, "properties");
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;

    #[test]
    fn manifest_declares_what_the_code_reads() {
        let m: Value = serde_json::from_str(include_str!("../operator.json")).unwrap();
        let props: Vec<&Value> = m["properties"].as_array().unwrap().iter().collect();
        let by_name = |n: &str| {
            props
                .iter()
                .find(|p| p["name"] == n)
                .unwrap_or_else(|| panic!("property '{n}' is not declared in operator.json"))
        };
        let d = Settings::default();
        assert_eq!(
            by_name("selectionQuantile")["defaultValue"],
            d.selection_quantile
        );
        assert_eq!(
            by_name("depthScoreThreshold")["defaultValue"],
            d.depth_score_threshold
        );
        assert_eq!(by_name("seed")["defaultValue"], d.seed as f64);
        assert_eq!(by_name("sample_factor")["defaultValue"], "");
        assert_eq!(by_name("threads")["defaultValue"], 0.0);
        assert_eq!(by_name("histograms")["defaultValue"], d.histograms);
        let allowed = [
            "OperatorSpec",
            "CrosstabSpec",
            "MetaFactor",
            "AxisSpec",
            "OperatorJoinSpec",
            "JoinOperator",
            "ColumnPair",
            "TableRelation",
            "Attribute",
            "Pair",
        ];
        fn walk(v: &Value, allowed: &[&str]) {
            match v {
                Value::Object(o) => {
                    if let Some(Value::String(k)) = o.get("kind") {
                        assert!(allowed.contains(&k.as_str()), "invented kind '{k}'");
                    }
                    o.values().for_each(|c| walk(c, allowed));
                }
                Value::Array(a) => a.iter().for_each(|c| walk(c, allowed)),
                _ => {}
            }
        }
        walk(&m["operatorSpec"], &allowed);
        let c = m["container"].as_str().unwrap();
        let tag = c.rsplit_once(':').map(|(_, t)| t).unwrap_or("");
        assert!(
            tag.split('.').count() == 3
                && tag
                    .split('.')
                    .all(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit())),
            "container must pin an exact x.y.z version tag"
        );
        assert_eq!(m["communicationProtocol"], "grpc");
    }
}
