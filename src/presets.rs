use anyhow::{Context, Result};
use clap::ValueEnum;

use crate::model::LearnedModel;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum ModelPreset {
    /// Built-in E. coli model with Dam (GATC) and Dcm (CCWGG) methylation
    #[value(name = "ecoli")]
    EColiG6mATCC5mCWGG,
}

impl ModelPreset {
    pub fn label(&self) -> &'static str {
        match self {
            ModelPreset::EColiG6mATCC5mCWGG => "E. coli (G6mATC + C5mCWGG)",
        }
    }

    pub fn cli_token(&self) -> &'static str {
        match self {
            ModelPreset::EColiG6mATCC5mCWGG => "ecoli",
        }
    }

    pub fn load_model(&self) -> Result<LearnedModel> {
        match self {
            ModelPreset::EColiG6mATCC5mCWGG => {
                let raw = include_str!("../models/e_coli_G6mATC_C5mCWGG.toml");
                serde_json::from_str(raw).context("Failed to parse bundled E. coli model")
            }
        }
    }
}
