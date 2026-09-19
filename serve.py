# serve.py — FastAPI service for the deployment-selected clinical Cox model

import json
import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from typing import Literal

# ---- Load exported model artifacts ----
COEF  = json.load(open("models/coef.json"))
INFO  = json.load(open("models/model_info.json"))
BASE  = pd.read_csv("models/baseline_surv.csv").sort_values("time").reset_index(drop=True)

MEAN_AGE_TRAIN = INFO["mean_age_train"]
GAMMA          = INFO["gamma"]  # =1, documented in 12_export_model.R as correct-by-construction

# Linear predictor AT the reference row used to build baseline_surv.csv
# (mean age, gender=female -> gendermale=0, stage_grouped=I -> all dummies=0)
LP_REFERENCE = COEF["age_years"] * MEAN_AGE_TRAIN

app = FastAPI(
    title="TCGA-LUAD Clinical Survival Model",
    description="Research/portfolio use only. Not for clinical decision-making. "
                "See /model-info for known limitations disclosed at export time."
)

class PatientInput(BaseModel):
    age_years: float = Field(..., ge=0, le=120, description="Age in years")
    gender: Literal["female", "male"]
    stage_grouped: Literal["I", "II", "III", "IV"] = Field(
        ..., description="Collapsed AJCC stage (I/II/III/IV) — see model_info.json for why 8-level staging was rejected")

class RiskOutput(BaseModel):
    linear_predictor: float
    risk_relative_to_reference: str
    survival_1yr: float
    survival_3yr: float
    survival_5yr: float
    held_out_test_cindex: float
    known_limitations: list
    disclaimer: str

def compute_lp(age_years: float, gender: str, stage_grouped: str) -> float:
    """Linear predictor using EXACTLY the 5 exported coefficients — one-hot vs 'I'/'female' reference."""
    lp = COEF["age_years"] * age_years
    lp += COEF["gendermale"] * (1 if gender == "male" else 0)
    lp += COEF.get("stage_groupedII",  0) * (1 if stage_grouped == "II"  else 0)
    lp += COEF.get("stage_groupedIII", 0) * (1 if stage_grouped == "III" else 0)
    lp += COEF.get("stage_groupedIV",  0) * (1 if stage_grouped == "IV"  else 0)
    return lp

def s_baseline(t: float) -> float:
    """Step-function lookup: S_baseline(t) at the reference patient, as computed by survfit() in R."""
    s = BASE.loc[BASE["time"] <= t, "S_mean"]
    return float(s.iloc[-1]) if len(s) else 1.0

def survival_at(t: float, lp_patient: float) -> float:
    """S(t | patient) = S_baseline(t) ^ exp(GAMMA * (lp_patient - lp_reference))"""
    exponent = np.exp(GAMMA * (lp_patient - LP_REFERENCE))
    return float(s_baseline(t) ** exponent)

@app.post("/predict", response_model=RiskOutput)
def predict(p: PatientInput):
    lp_patient = compute_lp(p.age_years, p.gender, p.stage_grouped)
    rel = lp_patient - LP_REFERENCE
    risk_label = (
        f"{'Higher' if rel > 0 else 'Lower' if rel < 0 else 'Equal'} risk than reference "
        f"(mean age {MEAN_AGE_TRAIN:.1f}, female, Stage I) — "
        f"hazard ratio vs reference = {np.exp(rel):.2f}"
    )
    return RiskOutput(
        linear_predictor=lp_patient,
        risk_relative_to_reference=risk_label,
        survival_1yr=survival_at(365, lp_patient),
        survival_3yr=survival_at(1095, lp_patient),
        survival_5yr=survival_at(1825, lp_patient),
        held_out_test_cindex=INFO["held_out_test_cindex"],
        known_limitations=INFO["known_limitations"],
        disclaimer="Research use only. Not for clinical decision-making. "
                   "This model prioritizes statistical validity (proportional-hazards "
                   "compliance) over raw discrimination — see /model-info for the full "
                   "model-selection rationale.",
    )

@app.get("/model-info")
def model_info():
    """Full transparency endpoint — surfaces every disclosed limitation, not just a summary."""
    return INFO

@app.get("/health")
def health():
    return {"status": "ok", "model_type": INFO["model_type"], "n_coefficients": len(COEF)}