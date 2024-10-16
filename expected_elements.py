covalent_radius = {  # from wikipedia
    1: 0.31,
    5: 0.84,
    6: 0.76,
    7: 0.71,
    8: 0.66,
    9: 0.57,
    12: 0.00,  # hack to avoid bonds with metals
    14: 1.11,
    15: 1.07,
    16: 1.05,
    17: 1.02,
    # 19: 2.03,
    20: 0.00,
    # 24: 1.39,
    25: 0.00,  # hack to avoid bonds with metals
    26: 0.00,
    30: 0.00,  # hack to avoid bonds with metals
    # 34: 1.20,
    35: 1.20,
    53: 1.39,
}
list_of_AD_elements_as_AtomicNum = list(covalent_radius.keys())