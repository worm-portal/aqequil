from .AqSpeciation import AqEquil, load, compare, Speciation
from .MassTransfer import Mass_Transfer, react, Reactant, Gas, Prepare_Reaction, Mixing_Fluid, join_mixes


def _configure_plotly_renderer():
    """
    Auto-detect the notebook environment and configure Plotly's renderer.
    This ensures plots display correctly in JupyterLab, classic notebooks,
    VS Code, Google Colab, and other environments.
    """
    try:
        import plotly.io as pio
        import os
        import sys

        # Check if we're in an IPython/Jupyter environment
        try:
            from IPython import get_ipython
            ipython = get_ipython()
            if ipython is None:
                # Not in IPython, leave default (browser)
                return
        except ImportError:
            return

        # Google Colab
        if 'google.colab' in sys.modules:
            pio.renderers.default = 'colab'
            return

        # Kaggle
        if 'KAGGLE_KERNEL_RUN_TYPE' in os.environ:
            pio.renderers.default = 'kaggle'
            return

        # VS Code notebooks
        if 'VSCODE_PID' in os.environ or 'VSCODE_CWD' in os.environ:
            pio.renderers.default = 'notebook'
            return

        # Databricks
        if 'DATABRICKS_RUNTIME_VERSION' in os.environ:
            pio.renderers.default = 'databricks'
            return

        # Check if running in a Jupyter kernel
        if ipython is not None and hasattr(ipython, 'kernel'):
            # Try to detect JupyterLab vs classic notebook
            # JupyterLab sets JUPYTERHUB_* or has specific config
            if 'JUPYTERHUB_API_TOKEN' in os.environ:
                # JupyterHub (often JupyterLab)
                pio.renderers.default = 'iframe'
                return

            # Check for JupyterLab by looking at the frontend
            try:
                # In JupyterLab, this environment variable may be set
                if 'JPY_SESSION_NAME' in os.environ:
                    # Could be either - iframe works for both
                    pio.renderers.default = 'iframe'
                    return
            except:
                pass

            # Default for Jupyter environments: iframe is most reliable
            pio.renderers.default = 'iframe'
            return

    except Exception:
        # If anything fails, silently continue with Plotly's default
        pass


# Configure Plotly renderer on import
_configure_plotly_renderer()