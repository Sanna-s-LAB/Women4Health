# Python script
# This script contains the object used in the script n.4 of this pipeline
# to run compute the residuals
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 13/01/2025
# Python version: v3.10.12






import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.metrics import r2_score, root_mean_squared_error
from sklearn.inspection import permutation_importance
from sklearn.base import clone




class GeneralizedModel:
    def __init__(self, data, outcomes, covariates, model, param_grid =None):
        """
        Initialize the GeneralizedModel class.

        Parameters:
        - data: pandas DataFrame containing the dataset.
        - outcomes: list of outcome variables.
        - covariates: list of covariate variables.
        """
        self.data = data
        self.outcomes = outcomes
        self.covariates = covariates
        self.formula = {out : out + ' ~ ' + ' + '.join(covariates) for out in outcomes}
        self.model = model
        self.param_grid = param_grid
        self.models = None
        self.residuals = None
        self.r_squared = None
        self.error = None

    def _prepare_data(self, specific_outcome):
        """
        Prepare the dataset for model training and testing.
        """
        df = self.data[self.outcomes + self.covariates]
        df = df[~df[specific_outcome].isna()]  # Drop rows with missing outcome
        X = df[self.covariates]
        y = df[specific_outcome]

        # Handle categorical variables
        categorical_columns = X.select_dtypes(include=['object', 'category']).columns
        X = pd.get_dummies(X, columns=categorical_columns, drop_first=True)

        return X, y


    def fit(self):
        """
        Fit models to predict outcomes using the specified covariates.
        """
        r_squared = {}
        err = {}
        impo = {}

        for specific_outcome in self.outcomes:

            X, y = self._prepare_data(specific_outcome)
            # Split data
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.10, random_state=42)

            # Train model without cross-validation if param_grid is None
            if self.param_grid:
                grid_search = GridSearchCV(
                    estimator=self.model,
                    param_grid=self.param_grid,
                    scoring='neg_mean_squared_error',
                    cv=10,
                    n_jobs=-1
                )
                grid_search.fit(X_train, y_train)
                best_model = grid_search.best_estimator_
            else:
                best_model = self.model
                best_model.fit(X_train, y_train)

            pred = best_model.predict(X_test)
            r_squared[specific_outcome] = r2_score(y_test, pred)
            err[specific_outcome] = root_mean_squared_error(y_test, pred)
            
            # Feature importance
            result = permutation_importance(best_model, X, y, n_repeats=100, random_state=42, n_jobs=-1)
            impo[specific_outcome] = result.importances_mean

        self.r_squared = pd.DataFrame.from_dict(r_squared, orient='index', columns=['R²_score']).reset_index()
        self.r_squared.columns = ['Feature', 'R²_score']

        self.error = pd.DataFrame.from_dict(err, orient='index', columns=['Mean_squared_error']).reset_index()
        self.error.columns = ['Feature', 'Mean_squared_error']

        self.importance = pd.DataFrame.from_dict(impo, orient='index', columns=X.columns).reset_index()
        self.importance = pd.melt(self.importance, id_vars=['index'], var_name='Variable', value_name='Importance')
        self.importance.columns = ['Feature', 'Variable', 'Importance']

    def compute_resid(self):
        residuals_df = []
        mdl = {}

        for specific_outcome in self.outcomes:
            # Prepare data
            X, y = self._prepare_data(specific_outcome)
            base_model = clone(self.model)
            # Fit model
            if self.param_grid:
                grid_search = GridSearchCV(
                    estimator=base_model,
                    param_grid=self.param_grid,
                    scoring='neg_mean_squared_error',
                    cv=10,
                    n_jobs=-1,
                )
                grid_search.fit(X, y)
                best_model = grid_search.best_estimator_
            else:
                best_model = base_model
                best_model.fit(X, y)

            # Compute residuals
            residuals_values = y - best_model.predict(X)
            temp_df = pd.DataFrame({
                'Outcome': specific_outcome,
                'Residual': residuals_values.values,
                'Observation': y.index
            })
            residuals_df.append(temp_df)

            # Save the model using specific_outcome as the key
            mdl[specific_outcome] = best_model

        # Save models and residuals
        self.models = mdl
        self.residuals = pd.concat(residuals_df, ignore_index=True).sort_values(by=['Outcome', 'Observation'])


    def compute_p_values(self, n_permutations=100):
        """
        Compute permutation-based p-values for feature importance.
        """
        if not hasattr(self, 'importance'):
            print("Model must be fitted before computing p-values.")
            return

        p_values = {}

        for specific_outcome in self.outcomes:
            X, y = self._prepare_data(specific_outcome)

            # Reference model
            ref_model = self.models[specific_outcome]
            original_r2 = r2_score(y, ref_model.predict(X))

            p_values_outcome = {}
            for cov in X.columns:
                permuted_r2_scores = []

                for _ in range(n_permutations):
                    permuted_X = X.copy()
                    permuted_X[cov] = np.random.permutation(permuted_X[cov])
                    permuted_r2 = r2_score(y, ref_model.predict(permuted_X))
                    permuted_r2_scores.append(permuted_r2)

                # Calculate p-value
                p_value = np.mean(np.array(permuted_r2_scores) >= original_r2)
                p_values_outcome[cov] = p_value

            p_values[specific_outcome] = p_values_outcome

        # Create a DataFrame for p-values
        self.p_values = pd.DataFrame.from_dict(p_values, orient='index').reset_index()
        self.p_values = pd.melt(self.p_values, id_vars=['index'], var_name='Variable', value_name='p_value')
        self.p_values.columns = ['Feature', 'Variable', 'p_value']

    def print_res(self, res):

        if hasattr(self, res) and getattr(self, res) is not None:     
            if res == 'importance':
                plt.figure(figsize=(12, 6))
                sns.barplot(data=self.importance, hue='Variable', y='Importance', x='Feature', errorbar=('ci', False))
                plt.xticks(rotation=90, fontsize=10)
                plt.xlabel("Feature")
                plt.ylabel("Importance")
                plt.title("Feature Importances by Outcome")
                plt.legend(bbox_to_anchor=(1.2, 1), borderaxespad=0)
                plt.tight_layout()
                plt.show()
                    
            elif res == 'r_squared':
                plt.figure(figsize=(8, 5))
                sns.barplot(data=self.r_squared, x='Feature', y='R²_score', errorbar=('ci', False))
                plt.xticks(rotation=90, fontsize=10)
                plt.xlabel("Outcome")
                plt.ylabel("R² Score")
                plt.title("R² Scores by Outcome")
                plt.tight_layout()
                plt.show()      
                
            elif res == 'p_values':
                plt.figure(figsize=(12, 6))
                ax = sns.barplot(data=self.p_values, x='Variable', y='p_value', hue = 'Feature', errorbar=('ci',False))
                # Define a significance threshold for testing
                significance_threshold = 0.05  # For debugging, set this to 0.5 to see if asterisks appear
                
                plt.axhline(y=significance_threshold, color='red', linestyle='--', linewidth=1, label=f"Significance Threshold (p={significance_threshold})")
                plt.xticks(rotation=90, fontsize=10)
                plt.xlabel("Feature")
                plt.ylabel("p-value")
                plt.title("P-values by Feature and Hormone")
                plt.legend()
                
        else:
            print('You have to fit the model before printing or the requested result is unavailable.')

    def plot_density(self, columns, separate=False):
        """
        Plot the density of specified numerical columns in the dataset.

        Parameters:
        - columns: list of column names to plot.
        - separate: bool, if True, creates separate subplots for each column.
        """
        # Check if all columns exist in the dataset
        missing_columns = [col for col in columns if col not in self.data.columns]
        if missing_columns:
            print(f"The following columns are not in the dataset: {missing_columns}")
            return
        
        # Filter for numerical columns
        numerical_columns = [col for col in columns if self.data[col].dtype in ['float64', 'int64']]
        if not numerical_columns:
            print("No numerical columns found in the provided list.")
            return

        if separate:
            # Create subplots for each column
            n_cols = 2  # Number of columns in the subplot grid
            n_rows = -(-len(numerical_columns) // n_cols)  # Calculate rows needed
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 6 * n_rows), constrained_layout=True)
            axes = axes.flatten()

            for idx, column in enumerate(numerical_columns):
                sns.kdeplot(self.data[column], ax=axes[idx], fill=True, color='blue')
                axes[idx].set_title(f"Density Plot: {column}")
                axes[idx].set_xlabel("Value")
                axes[idx].set_ylabel("Density")
            
            # Hide any unused subplots
            for ax in axes[len(numerical_columns):]:
                ax.set_visible(False)
            
            plt.show()

        else:
            # Single plot for all densities
            plt.figure(figsize=(10, 6))
            for column in numerical_columns:
                sns.kdeplot(self.data[column], label=column, fill=True)

            plt.title("Density Plot of Selected Columns")
            plt.xlabel("Value")
            plt.ylabel("Density")
            plt.legend(title="Columns")
            plt.tight_layout()
            plt.show()
