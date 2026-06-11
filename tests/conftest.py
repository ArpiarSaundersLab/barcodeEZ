import pytest
import pandas as pd
from pathlib import Path

REFERENCE_DATA = Path(__file__).parent / 'data'


@pytest.fixture(scope='session')
def sdf2_ref():
    return pd.read_csv(REFERENCE_DATA / 'SDF2_260331_IDTOrder_Sites1-3_Compiled_v2.csv')


@pytest.fixture(scope='session')
def sdf3_ref():
    return pd.read_csv(REFERENCE_DATA / 'SDF3_260406_IDTOrder_Sites1-4_Compiled_v2.csv')
