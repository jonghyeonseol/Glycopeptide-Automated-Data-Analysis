"""
Statistical Validation Module

Provides rigorous statistical validation methods for glycoproteomics analysis,
including bootstrap VIP validation, cross-validation, effect size calculations,
and permutation tests for publication-quality results.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-06
"""

import numpy as np
import pandas as pd
from typing import List
from dataclasses import dataclass
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report
import warnings


@dataclass
class BootstrapVIPResult:
    """Results from bootstrap VIP validation"""
    vip_mean: pd.Series
    vip_std: pd.Series
    vip_ci_lower: pd.Series
    vip_ci_upper: pd.Series
    stability_score: pd.Series  # Proportion of times VIP > 1.0
    n_iterations: int
    feature_names: List[str]

    def get_stable_biomarkers(self, stability_threshold: float = 0.8,
                              vip_threshold: float = 1.0) -> pd.DataFrame:
        """
        Identify stable biomarkers based on bootstrap results

        Parameters:
        -----------
        stability_threshold : float
            Minimum stability score (default: 0.8 = 80% of iterations)
        vip_threshold : float
            Minimum VIP threshold (default: 1.0)

        Returns:
        --------
        pd.DataFrame with stable biomarkers and their statistics
        """
        stable_mask = (self.stability_score >= stability_threshold) & \
            (self.vip_mean >= vip_threshold)

        results = pd.DataFrame({
            'Feature': self.feature_names,
            'VIP_Mean': self.vip_mean,
            'VIP_Std': self.vip_std,
            'VIP_CI_Lower': self.vip_ci_lower,
            'VIP_CI_Upper': self.vip_ci_upper,
            'Stability_Score': self.stability_score
        })

        return results[stable_mask].sort_values('VIP_Mean', ascending=False)


@dataclass
class CrossValidationResult:
    """Results from PLS-DA cross-validation"""
    accuracy_scores: np.ndarray
    accuracy_mean: float
    accuracy_std: float
    roc_auc_scores: np.ndarray
    roc_auc_mean: float
    roc_auc_std: float
    n_folds: int
    n_components: int
    classification_reports: List[str]

    def summary(self) -> str:
        """Generate summary report"""
        report = []
        report.append(f"Cross-Validation Results ({self.n_folds}-Fold)")
        report.append(f"PLS-DA Components: {self.n_components}")
        report.append(f"\nAccuracy: {self.accuracy_mean:.3f} ± {self.accuracy_std:.3f}")
        report.append(f"ROC-AUC: {self.roc_auc_mean:.3f} ± {self.roc_auc_std:.3f}")
        report.append(f"\nFold-wise Accuracy: {', '.join([f'{x:.3f}' for x in self.accuracy_scores])}")
        report.append(f"Fold-wise ROC-AUC: {', '.join([f'{x:.3f}' for x in self.roc_auc_scores])}")
        return '\n'.join(report)


@dataclass
class EffectSizeResult:
    """Results from effect size calculations"""
    cohens_d: pd.Series
    effect_magnitude: pd.Series  # small/medium/large
    group1_mean: pd.Series
    group2_mean: pd.Series
    group1_std: pd.Series
    group2_std: pd.Series
    pooled_std: pd.Series

    def get_large_effects(self, threshold: float = 0.8) -> pd.DataFrame:
        """Get features with large effect sizes (|d| >= threshold)"""
        large_mask = np.abs(self.cohens_d) >= threshold

        results = pd.DataFrame({
            'Feature': self.cohens_d.index,
            'Cohens_d': self.cohens_d,
            'Effect_Magnitude': self.effect_magnitude,
            'Group1_Mean': self.group1_mean,
            'Group2_Mean': self.group2_mean,
            'Pooled_Std': self.pooled_std
        })

        return results[large_mask].sort_values('Cohens_d',
                                               key=lambda x: np.abs(x),
                                               ascending=False)


@dataclass
class PermutationTestResult:
    """Results from permutation test"""
    observed_statistic: float
    null_distribution: np.ndarray
    p_value: float
    n_permutations: int
    statistic_name: str

    def is_significant(self, alpha: float = 0.05) -> bool:
        """Check if result is significant at given alpha level"""
        return self.p_value < alpha


class StatisticalValidator:
    """
    Statistical validation for glycoproteomics analysis

    Provides publication-quality statistical validation including:
    - Bootstrap VIP validation with confidence intervals
    - PLS-DA cross-validation with performance metrics
    - Cohen's d effect size calculations
    - Permutation tests for PCA/PLS-DA significance
    """

    def __init__(self, random_state: int = 42):
        """
        Initialize validator

        Parameters:
        -----------
        random_state : int
            Random seed for reproducibility
        """
        self.random_state = random_state
        np.random.seed(random_state)

    def bootstrap_vip_validation(self,
                                 X: np.ndarray,
                                 y: np.ndarray,
                                 feature_names: List[str],
                                 n_iterations: int = 1000,
                                 n_components: int = 2,
                                 confidence_level: float = 0.95) -> BootstrapVIPResult:
        """
        Perform bootstrap validation of VIP scores

        Parameters:
        -----------
        X : np.ndarray
            Feature matrix (samples x features)
        y : np.ndarray
            Group labels (0/1)
        feature_names : List[str]
            Names of features
        n_iterations : int
            Number of bootstrap iterations (default: 1000)
        n_components : int
            Number of PLS components (default: 2)
        confidence_level : float
            Confidence level for intervals (default: 0.95)

        Returns:
        --------
        BootstrapVIPResult with mean, std, CI, and stability scores
        """
        n_samples, n_features = X.shape
        vip_scores_all = np.zeros((n_iterations, n_features))

        print("\n🔬 Bootstrap VIP Validation")
        print(f"   Iterations: {n_iterations}")
        print(f"   Samples: {n_samples}, Features: {n_features}")
        print(f"   Confidence Level: {confidence_level * 100:.0f}%")

        for i in range(n_iterations):
            if (i + 1) % 100 == 0:
                print(f"   Progress: {i + 1}/{n_iterations} iterations...", end='\r')

            # Bootstrap resampling with replacement
            indices = np.random.choice(n_samples, size=n_samples, replace=True)
            X_boot = X[indices, :]
            y_boot = y[indices]

            # Ensure both classes present in bootstrap sample
            if len(np.unique(y_boot)) < 2:
                # Resample until we get both classes
                while len(np.unique(y_boot)) < 2:
                    indices = np.random.choice(n_samples, size=n_samples, replace=True)
                    X_boot = X[indices, :]
                    y_boot = y[indices]

            # Calculate VIP scores for this iteration
            vip_scores = self._calculate_vip_scores(X_boot, y_boot, n_components)
            vip_scores_all[i, :] = vip_scores

        print(f"   ✓ Completed {n_iterations} iterations")

        # Calculate statistics
        vip_mean = np.mean(vip_scores_all, axis=0)
        vip_std = np.std(vip_scores_all, axis=0)

        # Confidence intervals
        alpha = 1 - confidence_level
        vip_ci_lower = np.percentile(vip_scores_all, alpha / 2 * 100, axis=0)
        vip_ci_upper = np.percentile(vip_scores_all, (1 - alpha / 2) * 100, axis=0)

        # Stability score: proportion of iterations where VIP > 1.0
        stability_score = np.mean(vip_scores_all > 1.0, axis=0)

        return BootstrapVIPResult(
            vip_mean=pd.Series(vip_mean, index=feature_names),
            vip_std=pd.Series(vip_std, index=feature_names),
            vip_ci_lower=pd.Series(vip_ci_lower, index=feature_names),
            vip_ci_upper=pd.Series(vip_ci_upper, index=feature_names),
            stability_score=pd.Series(stability_score, index=feature_names),
            n_iterations=n_iterations,
            feature_names=feature_names
        )

    def _calculate_vip_scores(self, X: np.ndarray, y: np.ndarray,
                              n_components: int = 2) -> np.ndarray:
        """Calculate VIP scores for PLS-DA model"""
        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Fit PLS model
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pls = PLSRegression(n_components=n_components, scale=False)
            pls.fit(X_scaled, y)

        # Calculate VIP scores
        t = pls.x_scores_
        w = pls.x_weights_
        q = pls.y_loadings_

        p, h = w.shape
        vips = np.zeros((p,))

        s = np.diag(np.dot(np.dot(np.dot(t.T, t), q.T), q)).reshape(h, -1)
        total_s = np.sum(s)

        for i in range(p):
            weight = np.array([(w[i, j] / np.linalg.norm(w[:, j]))**2
                              for j in range(h)])
            vips[i] = np.sqrt(p * (np.dot(s.T, weight)) / total_s)

        return vips

    def cross_validate_plsda(self,
                             X: np.ndarray,
                             y: np.ndarray,
                             n_components: int = 2,
                             n_folds: int = 10) -> CrossValidationResult:
        """
        Perform cross-validation for PLS-DA model

        Parameters:
        -----------
        X : np.ndarray
            Feature matrix (samples x features)
        y : np.ndarray
            Group labels (0/1)
        n_components : int
            Number of PLS components
        n_folds : int
            Number of cross-validation folds

        Returns:
        --------
        CrossValidationResult with accuracy, ROC-AUC, and detailed reports
        """
        print("\n🔬 PLS-DA Cross-Validation")
        print(f"   Folds: {n_folds}")
        print(f"   Components: {n_components}")

        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Setup cross-validation
        cv = StratifiedKFold(n_splits=n_folds, shuffle=True,
                             random_state=self.random_state)

        accuracy_scores = []
        roc_auc_scores = []
        classification_reports = []

        for fold, (train_idx, test_idx) in enumerate(cv.split(X_scaled, y), 1):
            X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Fit PLS-DA
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pls = PLSRegression(n_components=n_components, scale=False)
                pls.fit(X_train, y_train)

            # Predict
            y_pred_proba = pls.predict(X_test).ravel()
            y_pred = (y_pred_proba >= 0.5).astype(int)

            # Metrics
            accuracy = accuracy_score(y_test, y_pred)
            roc_auc = roc_auc_score(y_test, y_pred_proba)

            accuracy_scores.append(accuracy)
            roc_auc_scores.append(roc_auc)

            # Classification report
            report = classification_report(y_test, y_pred,
                                           target_names=['Normal', 'Cancer'],
                                           output_dict=False)
            classification_reports.append(f"Fold {fold}:\n{report}")

            print(f"   Fold {fold}: Accuracy={accuracy:.3f}, ROC-AUC={roc_auc:.3f}")

        accuracy_scores = np.array(accuracy_scores)
        roc_auc_scores = np.array(roc_auc_scores)

        print(f"   ✓ Mean Accuracy: {np.mean(accuracy_scores):.3f} ± {np.std(accuracy_scores):.3f}")
        print(f"   ✓ Mean ROC-AUC: {np.mean(roc_auc_scores):.3f} ± {np.std(roc_auc_scores):.3f}")

        return CrossValidationResult(
            accuracy_scores=accuracy_scores,
            accuracy_mean=np.mean(accuracy_scores),
            accuracy_std=np.std(accuracy_scores),
            roc_auc_scores=roc_auc_scores,
            roc_auc_mean=np.mean(roc_auc_scores),
            roc_auc_std=np.std(roc_auc_scores),
            n_folds=n_folds,
            n_components=n_components,
            classification_reports=classification_reports
        )

    def calculate_cohens_d(self,
                           group1: np.ndarray,
                           group2: np.ndarray,
                           feature_names: List[str]) -> EffectSizeResult:
        """
        Calculate Cohen's d effect size for each feature

        Parameters:
        -----------
        group1 : np.ndarray
            Data for group 1 (samples x features)
        group2 : np.ndarray
            Data for group 2 (samples x features)
        feature_names : List[str]
            Names of features

        Returns:
        --------
        EffectSizeResult with Cohen's d and effect magnitude
        """
        print("\n🔬 Cohen's d Effect Size Calculation")
        print(f"   Group 1: {group1.shape[0]} samples")
        print(f"   Group 2: {group2.shape[0]} samples")
        print(f"   Features: {len(feature_names)}")

        # Calculate means
        mean1 = np.nanmean(group1, axis=0)
        mean2 = np.nanmean(group2, axis=0)

        # Calculate standard deviations
        std1 = np.nanstd(group1, axis=0, ddof=1)
        std2 = np.nanstd(group2, axis=0, ddof=1)

        # Calculate pooled standard deviation
        n1 = group1.shape[0]
        n2 = group2.shape[0]
        pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))

        # Calculate Cohen's d
        cohens_d = (mean1 - mean2) / pooled_std

        # Classify effect magnitude
        effect_magnitude = pd.Series(['negligible'] * len(cohens_d), index=feature_names)
        effect_magnitude[np.abs(cohens_d) >= 0.2] = 'small'
        effect_magnitude[np.abs(cohens_d) >= 0.5] = 'medium'
        effect_magnitude[np.abs(cohens_d) >= 0.8] = 'large'

        # Count effect sizes
        n_large = np.sum(np.abs(cohens_d) >= 0.8)
        n_medium = np.sum((np.abs(cohens_d) >= 0.5) & (np.abs(cohens_d) < 0.8))
        n_small = np.sum((np.abs(cohens_d) >= 0.2) & (np.abs(cohens_d) < 0.5))

        print(f"   ✓ Effect Sizes: {n_large} large, {n_medium} medium, {n_small} small")

        return EffectSizeResult(
            cohens_d=pd.Series(cohens_d, index=feature_names),
            effect_magnitude=effect_magnitude,
            group1_mean=pd.Series(mean1, index=feature_names),
            group2_mean=pd.Series(mean2, index=feature_names),
            group1_std=pd.Series(std1, index=feature_names),
            group2_std=pd.Series(std2, index=feature_names),
            pooled_std=pd.Series(pooled_std, index=feature_names)
        )

    def permutation_test_pca(self,
                             X: np.ndarray,
                             y: np.ndarray,
                             n_permutations: int = 1000) -> PermutationTestResult:
        """
        Permutation test for PCA group separation

        Tests whether observed between-group separation in PC1-PC2 space
        is greater than expected by chance.

        Parameters:
        -----------
        X : np.ndarray
            Feature matrix (samples x features)
        y : np.ndarray
            Group labels (0/1)
        n_permutations : int
            Number of permutations

        Returns:
        --------
        PermutationTestResult with p-value and null distribution
        """

        print("\n🔬 PCA Permutation Test")
        print(f"   Permutations: {n_permutations}")

        # Calculate observed statistic (between-group distance in PC space)
        observed_stat = self._calculate_pca_separation(X, y)
        print(f"   Observed separation: {observed_stat:.4f}")

        # Generate null distribution
        null_distribution = np.zeros(n_permutations)

        for i in range(n_permutations):
            if (i + 1) % 100 == 0:
                print(f"   Progress: {i + 1}/{n_permutations} permutations...", end='\r')

            # Permute labels
            y_perm = np.random.permutation(y)

            # Calculate statistic under null
            null_distribution[i] = self._calculate_pca_separation(X, y_perm)

        print(f"   ✓ Completed {n_permutations} permutations")

        # Calculate p-value
        p_value = np.mean(null_distribution >= observed_stat)

        print(f"   ✓ P-value: {p_value:.4f}")

        return PermutationTestResult(
            observed_statistic=observed_stat,
            null_distribution=null_distribution,
            p_value=p_value,
            n_permutations=n_permutations,
            statistic_name="PCA_separation"
        )

    def _calculate_pca_separation(self, X: np.ndarray, y: np.ndarray) -> float:
        """Calculate between-group separation in PCA space"""
        from sklearn.decomposition import PCA

        # Standardize and perform PCA
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_scaled)

        # Calculate centroids for each group
        group0_centroid = np.mean(X_pca[y == 0], axis=0)
        group1_centroid = np.mean(X_pca[y == 1], axis=0)

        # Euclidean distance between centroids
        separation = np.linalg.norm(group1_centroid - group0_centroid)

        return separation
