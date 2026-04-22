import pytest
import numpy as np
from src.pcr_autolearner_ai import PCRCluster

@pytest.fixture
def sample_data():
    """Create synthetic PCR data for testing clustering.
    All data shifted by +100 to ensure non-zero start for normalization test.
    """
    # Create 3 distinct clusters
    t = np.linspace(0, 45, 45)
    
    # Baseline offset
    baseline = 100.0
    
    # Cluster 1: Sigmoidal curve (Growth)
    # Range 100 to ~101 or more depending on scale. Let's make it more distinct
    # Sigmoid 0 to 1 -> * 100 -> 0 to 100. Total 100 to 200.
    c1 = 100 * (1 / (1 + np.exp(-0.5 * (t - 20))))
    data1 = np.tile(c1, (10, 1)) + np.random.normal(0, 1, (10, 45)) + baseline
    
    # Cluster 2: Flat line (No growth)
    c2 = np.zeros_like(t)
    data2 = np.tile(c2, (10, 1)) + np.random.normal(0, 1, (10, 45)) + baseline
    
    # Cluster 3: Linear growth (Background / Drift)
    c3 = 2.0 * t  # 0 to 90
    data3 = np.tile(c3, (10, 1)) + np.random.normal(0, 1, (10, 45)) + baseline
    
    return np.vstack([data1, data2, data3])

def test_kmeans_clustering(sample_data):
    """Test KMeans clustering."""
    clusterer = PCRCluster(n_clusters=3, method="kmeans")
    normalized_data, _ = clusterer.preprocess(sample_data)
    labels = clusterer.fit_predict(normalized_data)
    
    assert len(labels) == 30
    assert len(np.unique(labels)) == 3

def test_gmm_clustering(sample_data):
    """Test GMM clustering."""
    clusterer = PCRCluster(n_clusters=3, method="gmm")
    normalized_data, _ = clusterer.preprocess(sample_data)
    labels = clusterer.fit_predict(normalized_data)
    
    assert len(labels) == 30
    assert len(np.unique(labels)) == 3

def test_bgmm_clustering(sample_data):
    """Test BGMM clustering."""
    clusterer = PCRCluster(n_clusters=5, method="bgmm")
    normalized_data, _ = clusterer.preprocess(sample_data)
    labels = clusterer.fit_predict(normalized_data)
    
    assert len(labels) == 30
    n_labels = len(np.unique(labels))
    assert n_labels >= 2 and n_labels <= 5

def test_preprocess_filtering():
    """Test that preprocess filters out samples starting with 0."""
    clusterer = PCRCluster()
    
    # Data 1: starts with 0
    d1 = np.array([0, 10, 20])
    # Data 2: starts with 10
    d2 = np.array([10, 20, 30])
    
    data = np.vstack([d1, d2])
    
    normalized_data, valid_indices = clusterer.preprocess(data)
    
    assert len(normalized_data) == 1
    assert len(valid_indices) == 1
    assert valid_indices[0] == 1 # The second row (index 1) should remain
    
    # Check normalization: 10, 20, 30 divided by 10 -> 1, 2, 3
    expected = np.array([[1.0, 2.0, 3.0]])
    np.testing.assert_array_almost_equal(normalized_data, expected)

def test_preprocess_skip_normalization():
    """Test skipping normalization."""
    clusterer = PCRCluster()
    data = np.array([[0.0, 0.5, 1.0], [0.2, 0.4, 0.6]])
    
    normalized_data, valid_indices = clusterer.preprocess(data, skip_normalization=True)
    
    # Should be identical
    np.testing.assert_array_equal(normalized_data, data)
    assert len(valid_indices) == 2
    assert valid_indices.tolist() == [0, 1]

