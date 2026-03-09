import React, { useEffect, useState } from 'react';
import { Activity, Moon, Sun } from 'lucide-react';
import { MODEL_VERSION } from './constants';
import { SingleWorkflowPage } from './components/SingleWorkflowPage';
import { BatchWorkflowPage } from './components/BatchWorkflowPage';

const App: React.FC = () => {
  const [isDark, setIsDark] = useState<boolean>(true);
  const [mode, setMode] = useState<'single' | 'batch'>('single');

  useEffect(() => {
    if (isDark) {
      document.documentElement.classList.add('dark');
      document.body.classList.add('dark');
    } else {
      document.documentElement.classList.remove('dark');
      document.body.classList.remove('dark');
    }
  }, [isDark]);

  return (
    <div className="min-h-screen bg-slate-50 dark:bg-[#020617] transition-colors duration-300">
      <header className="sticky top-0 z-50 backdrop-blur-sm bg-white/80 dark:bg-slate-950/80 border-b border-slate-200 dark:border-slate-800">
        <div className="container mx-auto px-4 py-3 flex items-center justify-between gap-4">
          <div className="flex items-center gap-3">
            <div className="w-9 h-9 rounded-lg bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white shadow-lg shadow-cyan-500/30">
              <Activity className="w-5 h-5" />
            </div>
            <div>
              <div className="text-sm font-bold text-slate-900 dark:text-white tracking-tight">Virtual Cell AI</div>
              <div className="text-[11px] font-mono text-slate-500 dark:text-slate-400">{MODEL_VERSION}</div>
            </div>
          </div>

          <div className="flex items-center gap-1 rounded-xl border border-slate-200 dark:border-slate-800 bg-white/80 dark:bg-slate-900/80 p-1">
            <button
              type="button"
              onClick={() => setMode('single')}
              className={mode === 'single'
                ? 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all bg-cyan-50 dark:bg-cyan-950/40 text-cyan-700 dark:text-cyan-300 shadow-sm'
                : 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all text-slate-500 hover:text-slate-700 dark:text-slate-400 dark:hover:text-slate-200'}
            >
              Single Prediction
            </button>
            <button
              type="button"
              onClick={() => setMode('batch')}
              className={mode === 'batch'
                ? 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all bg-cyan-50 dark:bg-cyan-950/40 text-cyan-700 dark:text-cyan-300 shadow-sm'
                : 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all text-slate-500 hover:text-slate-700 dark:text-slate-400 dark:hover:text-slate-200'}
            >
              Batch Prediction
            </button>
          </div>

          <button
            onClick={() => setIsDark(!isDark)}
            className="p-2 rounded-lg border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-900 text-slate-700 dark:text-slate-300 hover:text-cyan-600 dark:hover:text-cyan-400 transition-colors"
            aria-label="Toggle theme"
          >
            {isDark ? <Sun className="w-4 h-4" /> : <Moon className="w-4 h-4" />}
          </button>
        </div>
      </header>

      <main className="container mx-auto px-4 py-8">
        {mode === 'single' ? (
          <SingleWorkflowPage isDark={isDark} />
        ) : (
          <BatchWorkflowPage isDark={isDark} />
        )}
      </main>
    </div>
  );
};

export default App;
