#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PET显像剂剂量计算器
用于根据患者体重计算推荐给药剂量
"""

# 显像剂剂量标准（MBq/kg 或固定剂量）
DOSE_STANDARDS = {
    "FDG": {
        "type": "weight_based",
        "adult_mbq_per_kg": 5.18,  # 0.15 mCi/kg ≈ 5.18 MBq/kg
        "adult_mci_per_kg": 0.15,
        "typical_range_mbq": (185, 370),  # 5-10 mCi
        "typical_range_mci": (5, 10),
        "pediatric_neuro_mbq": 96.2,  # 约2.6 mCi
        "pediatric_neuro_mci": 2.6,
        "notes": "成人按体重计算，儿童神经学检查用固定剂量"
    },
    "AV45": {
        "type": "fixed",
        "fixed_mbq": 370,  # 10 mCi
        "fixed_mci": 10,
        "notes": "固定剂量"
    },
    "AV1451": {
        "type": "fixed",
        "fixed_mbq": 370,  # 10 mCi
        "fixed_mci": 10,
        "notes": "固定剂量"
    },
    "FPCIT": {
        "type": "weight_or_fixed",
        "typical_mbq": 148,  # 3-5 mCi, 取中值4 mCi
        "typical_mci": 4,
        "range_mbq": (111, 185),
        "range_mci": (3, 5),
        "notes": "按体重或固定剂量；无需停药，辐射剂量极低"
    },
    "PSMA": {
        "type": "weight_based",
        "typical_mbq": 222,  # 4-8 mCi, 取中值6 mCi
        "typical_mci": 6,
        "range_mbq": (150, 300),
        "range_mci": (4, 8),
        "notes": "按体重计算"
    },
    "奥曲肽": {
        "type": "weight_based",
        "mbq_per_kg_range": (1.0, 2.7),  # 0.05-0.1 mCi/kg
        "mci_per_kg_range": (0.05, 0.1),
        "typical_range_mbq": (75, 200),
        "typical_range_mci": (2, 5),
        "notes": "按体重计算"
    },
    "FAPI": {
        "type": "fixed_or_weight",
        "typical_mbq": 148,  # 3-5 mCi, 取中值4 mCi
        "typical_mci": 4,
        "range_mbq": (111, 185),
        "range_mci": (3, 5),
        "notes": "固定或按体重"
    },
    "FES": {
        "type": "fixed",
        "fixed_mbq": 203,  # 5-6 mCi, 取中值5.5 mCi
        "fixed_mci": 5.5,
        "range_mbq": (185, 222),
        "range_mci": (5, 6),
        "notes": "固定剂量为主"
    },
    "DOPA": {
        "type": "weight_based",
        "mbq_per_kg": 3.7,  # 0.1 mCi/kg ≈ 3.7 MBq/kg
        "mci_per_kg": 0.1,
        "typical_range_mbq": (111, 185),
        "typical_range_mci": (3, 5),
        "notes": "按体重计算"
    }
}

# 核素半衰期（分钟）
HALF_LIVES = {
    "F-18": 109.8,  # 约110分钟
    "Ga-68": 68.0,  # 约68分钟
}

def calculate_dose(tracer, weight_kg=None, is_pediatric=False):
    """
    计算推荐给药剂量
    
    参数:
        tracer: 显像剂名称
        weight_kg: 患者体重(kg)，体重基础型显像剂需要
        is_pediatric: 是否为儿童
    
    返回:
        dict: 包含推荐剂量、范围、说明
    """
    tracer = tracer.upper()
    if tracer not in DOSE_STANDARDS:
        return {"error": f"不支持的显像剂: {tracer}"}
    
    standard = DOSE_STANDARDS[tracer]
    result = {
        "tracer": tracer,
        "type": standard["type"],
        "notes": standard["notes"]
    }
    
    if standard["type"] == "fixed":
        result["recommended_mbq"] = standard["fixed_mbq"]
        result["recommended_mci"] = standard["fixed_mci"]
        
    elif standard["type"] == "weight_based":
        if weight_kg is None:
            result["error"] = "该显像剂需要患者体重"
            result["typical_range_mbq"] = standard.get("typical_range_mbq")
            result["typical_range_mci"] = standard.get("typical_range_mci")
        else:
            if tracer == "FDG" and is_pediatric:
                result["recommended_mbq"] = standard["pediatric_neuro_mbq"]
                result["recommended_mci"] = standard["pediatric_neuro_mci"]
            elif "mbq_per_kg" in standard:
                mbq = standard["mbq_per_kg"] * weight_kg
                mci = standard["mci_per_kg"] * weight_kg
                result["recommended_mbq"] = round(mbq, 1)
                result["recommended_mci"] = round(mci, 2)
            else:
                # 使用范围中值
                result["recommended_mbq"] = standard.get("typical_mbq", standard["range_mbq"][1])
                result["recommended_mci"] = standard.get("typical_mci", standard["range_mci"][1])
    
    elif standard["type"] in ["weight_or_fixed", "fixed_or_weight"]:
        if weight_kg:
            # 简单按体重比例估算
            result["recommended_mbq"] = round(standard["typical_mbq"] * weight_kg / 70, 1)
            result["recommended_mci"] = round(standard["typical_mci"] * weight_kg / 70, 2)
        else:
            result["recommended_mbq"] = standard["typical_mbq"]
            result["recommended_mci"] = standard["typical_mci"]
    
    return result

def decay_correction(initial_mbq, elapsed_minutes, nuclide="F-18"):
    """
    放射性衰变校正
    
    参数:
        initial_mbq: 初始活度(MBq)
        elapsed_minutes: 经过时间(分钟)
        nuclide: 核素类型
    
    返回:
        float: 校正后活度(MBq)
    """
    import math
    
    if nuclide not in HALF_LIVES:
        return None
    
    t_half = HALF_LIVES[nuclide]
    decayed_mbq = initial_mbq * math.exp(-0.693 * elapsed_minutes / t_half)
    return round(decayed_mbq, 2)

def format_dose_result(result):
    """格式化剂量计算结果"""
    if "error" in result:
        return f"计算错误: {result['error']}"
    
    output = f"【{result['tracer']} 剂量推荐】\n"
    
    if "recommended_mbq" in result:
        output += f"推荐剂量: {result['recommended_mbq']} MBq ({result['recommended_mci']} mCi)\n"
    
    if "typical_range_mbq" in result and result["typical_range_mbq"]:
        r = result["typical_range_mbq"]
        r_mci = result.get("typical_range_mci")
        output += f"常规范围: {r[0]}-{r[1]} MBq"
        if r_mci:
            output += f" ({r_mci[0]}-{r_mci[1]} mCi)"
        output += "\n"
    
    output += f"说明: {result['notes']}\n"
    return output

# 示例用法
if __name__ == "__main__":
    # 示例1: FDG成人，体重70kg
    result = calculate_dose("FDG", weight_kg=70)
    print(format_dose_result(result))
    print()
    
    # 示例2: FDG儿童神经检查
    result = calculate_dose("FDG", is_pediatric=True)
    print(format_dose_result(result))
    print()
    
    # 示例3: AV45固定剂量
    result = calculate_dose("AV45")
    print(format_dose_result(result))
    print()
    
    # 示例4: 衰变校正
    decayed = decay_correction(370, 60, "F-18")
    print(f"370 MBq F-18经过60分钟后剩余: {decayed} MBq")
