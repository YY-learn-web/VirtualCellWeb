
import React, { useState, useEffect } from 'react';
import { 
  Activity, 
  Dna, 
  Microscope, 
  Clock, 
  FlaskConical, 
  Upload, 
  ChevronDown,
  Cpu,
  Zap,
  Moon,
  Sun,
  Layers,
  FileSpreadsheet,
  Download,
  CheckCircle2
} from 'lucide-react';
import { Button } from './components/Button';
import { Input, RangeSlider } from './components/Input';
import { ResultsChart, SweepLineChart } from './components/ResultsChart';
import { MetabolomicsPage } from './components/MetabolomicsPage';
import { EnrichmentPage } from './components/EnrichmentPage';
import { 
  AppState, 
  PredictionRequest, 
  GeneExpressionResult, 
  CellLineOption, 
  IC50Prediction,
  IC50PredictionRequest,
  AnalysisMode,
  BatchMethod,
  SweepVariable,
  SweepPoint,
  BatchSweepRequest,
  BatchIC50Request,
  BatchIC50Response,
  BatchUploadData,
  SweepIC50Request,
  SweepIC50Response
} from './types';
import { PRESET_CELL_LINES, MAX_TIME_HOURS, MAX_DOSE_UM, MODEL_VERSION } from './constants';
// 现在使用的是api2.ts 
import { ASCENDService } from './services/api2';

const App: React.FC = () => {
  // --- Theme State ---
  const [isDark, setIsDark] = useState<boolean>(true);
  const [activePage, setActivePage] = useState<'transcriptomics' | 'metabolomics' | 'enrichment'>('transcriptomics');

  // --- App Logic State ---
  const [mode, setMode] = useState<AnalysisMode>('single');
  const [appState, setAppState] = useState<AppState>(AppState.IDLE);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);

  // --- Single Mode Inputs ---
  const [ic50Result, setIc50Result] = useState<IC50Prediction | null>(null);
  const [ic50Loading, setIc50Loading] = useState(false);
  const [ic50Error, setIc50Error] = useState<string | null>(null);
  const [ic50Display, setIc50Display] = useState<'ic50' | 'lnIc50'>('ic50');
  const [batchIc50Result, setBatchIc50Result] = useState<BatchIC50Response['data'] | null>(null);
  const [batchIc50Loading, setBatchIc50Loading] = useState(false);
  const [batchIc50Error, setBatchIc50Error] = useState<string | null>(null);
  const [sweepIc50Result, setSweepIc50Result] = useState<SweepIC50Response['data'] | null>(null);
  const [sweepIc50Loading, setSweepIc50Loading] = useState(false);
  const [sweepIc50Error, setSweepIc50Error] = useState<string | null>(null);

  const handleIC50Prediction = async () => {
    if (!singleResults?.length) {
      setIc50Error('Please run Module 1 first to generate the expression profile.');
      return;
    }
    setIc50Loading(true);
    setIc50Error(null);
    setIc50Display('ic50');
    try {
      const payload: IC50PredictionRequest = {
        cellLineId: selectedCellLine,
        smiles,
        expression: singleResults.map(({ geneId, expressionLevel }) => ({
          geneId,
          expressionLevel,
        })),
        metadata: { timeHours: time, doseUm: dose },
      };
      const response = await ASCENDService.predictIC50(payload);
      setIc50Result(response.data);
    } catch (err: any) {
      setIc50Error(err.message ?? 'IC50 Prediction failure');
      setIc50Result(null);
    } finally {
      setIc50Loading(false);
    }
  };

  const handleBatchIC50Prediction = async () => {
    if (!uploadResult?.fileId) {
      setBatchIc50Error('Please complete batch processing first.');
      return;
    }
    setBatchIc50Loading(true);
    setBatchIc50Error(null);
    try {
      const payload: BatchIC50Request = {
        fileId: uploadResult.fileId,
        cellLineId: batchCellLine,
        previewCount: 5,
      };
      const response = await ASCENDService.predictBatchIC50(payload);
      setBatchIc50Result(response.data);
    } catch (err: any) {
      setBatchIc50Error(err.message ?? 'Batch IC50 prediction failed');
      setBatchIc50Result(null);
    } finally {
      setBatchIc50Loading(false);
    }
  };

  const handleSweepIC50Prediction = async () => {
    if (!sweepResults?.length) {
      setSweepIc50Error('Please run the sweep first to generate transcriptomic profiles.');
      return;
    }
    setSweepIc50Loading(true);
    setSweepIc50Error(null);
    try {
      const payload: SweepIC50Request = {
        smiles: batchSmiles,
        cellLineId: batchCellLine,
        sweepVariable: sweepVar,
        fixedParamValue: batchFixedParam,
        range: { start: sweepStart, end: sweepEnd, steps: sweepSteps },
        previewCount: 5,
      };
      const response = await ASCENDService.predictSweepIC50(payload);
      setSweepIc50Result(response.data);
    } catch (err: any) {
      setSweepIc50Error(err.message ?? 'Sweep IC50 prediction failed');
      setSweepIc50Result(null);
    } finally {
      setSweepIc50Loading(false);
    }
  };


  const [smiles, setSmiles] = useState<string>('CC(=O)OC1=CC=CC=C1C(=O)O');
  const [time, setTime] = useState<number>(24);
  const [dose, setDose] = useState<number>(10);
  const [cellLineMode, setCellLineMode] = useState<'preset' | 'custom'>('preset');
  const [selectedCellLine, setSelectedCellLine] = useState<string>(PRESET_CELL_LINES[0].id);
  const [singleResults, setSingleResults] = useState<GeneExpressionResult[] | null>(null);
  const [metabolomicsSeedFile, setMetabolomicsSeedFile] = useState<File | null>(null);
  const [metabolomicsSeedLabel, setMetabolomicsSeedLabel] = useState<string | null>(null);
  const [enrichmentSeedResults, setEnrichmentSeedResults] = useState<GeneExpressionResult[] | null>(null);
  const [enrichmentSeedLabel, setEnrichmentSeedLabel] = useState<string | null>(null);
  const [enrichmentLoading, setEnrichmentLoading] = useState<boolean>(false);
  const enrichmentContextLabel = enrichmentSeedLabel
    ? enrichmentSeedLabel
    : singleResults
      ? selectedCellLine + ' | ' + time + 'h | ' + dose + 'uM'
      : undefined;

  // --- Batch Mode Inputs ---
  const [batchMethod, setBatchMethod] = useState<BatchMethod>('sweep');
  const [sweepVar, setSweepVar] = useState<SweepVariable>('dose');
  // Sweep Ranges
  const [sweepStart, setSweepStart] = useState<number>(0);
  const [sweepEnd, setSweepEnd] = useState<number>(100);
  const [sweepSteps, setSweepSteps] = useState<number>(6);
  // Fixed params for sweep
  const [batchSmiles, setBatchSmiles] = useState<string>('CC(=O)OC1=CC=CC=C1C(=O)O');
  const [batchFixedParam, setBatchFixedParam] = useState<number>(24); // Time if Dose sweep, Dose if Time sweep
  const [batchCellLine, setBatchCellLine] = useState<string>(PRESET_CELL_LINES[0].id);
  
  const [sweepResults, setSweepResults] = useState<SweepPoint[] | null>(null);
  const [uploadResult, setUploadResult] = useState<BatchUploadData | null>(null);
// 修改：根据是否有文件显示不同界面（没有文件显示上传界面，有文件显示结果界面），上传前检查是否选择了文件，没有选择文件不允许上传，将文件对象传递给API服务，跟踪用户的文件选择状�?
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
//修改：增�?75-390�?button  528-545�?
/*删除内容<div className="border-2 border-dashed border-slate-200 dark:border-slate-700 rounded-xl p-8 text-center hover:border-cyan-400 dark:hover:border-cyan-500 transition-colors bg-white dark:bg-slate-950/50 cursor-pointer group">
<button type="button" className="text-cyan-600 hover:underline flex items-center gap-1 text-xs font-bold mt-2">
<Button className="flex items-center gap-2">*/

  // Toggle Theme Class on Body
  useEffect(() => {
    if (isDark) {
      document.documentElement.classList.add('dark');
    } else {
      document.documentElement.classList.remove('dark');
    }
  }, [isDark]);

  useEffect(() => {
    const syncFromHash = () => {
      const hash = window.location.hash.replace('#', '');
      const resolved =
        hash === 'metabolomics'
          ? 'metabolomics'
          : hash === 'enrichment'
            ? 'enrichment'
            : 'transcriptomics';
      setActivePage(resolved);
    };
    syncFromHash();
    window.addEventListener('hashchange', syncFromHash);
    return () => window.removeEventListener('hashchange', syncFromHash);
  }, []);

  // --- Handlers ---

  const handlePredict = async (e: React.FormEvent) => {
    e.preventDefault();
    setAppState(AppState.LOADING);
    setErrorMsg(null);
    setSingleResults(null);
    setEnrichmentSeedResults(null);
    setEnrichmentSeedLabel(null);
    setIc50Result(null);
    setIc50Error(null);
    setIc50Loading(false);
    setSweepResults(null);
    setUploadResult(null);
    setBatchIc50Result(null);
    setBatchIc50Error(null);
    setBatchIc50Loading(false);
    setSweepIc50Result(null);
    setSweepIc50Error(null);
    setSweepIc50Loading(false);

    try {
      if (mode === 'single') {
        const payload: PredictionRequest = {
          smiles,
          time,
          dose,
          cellLineMode,
          cellLineId: cellLineMode === 'preset' ? selectedCellLine : undefined,
        };
        const response = await ASCENDService.predictExpression(payload);
        setSingleResults(response.data);
      } else {
        // BATCH MODE
        if (batchMethod === 'sweep') {
           const payload: BatchSweepRequest = {
             smiles: batchSmiles,
             cellLineId: batchCellLine,
             sweepVariable: sweepVar,
             fixedParamValue: batchFixedParam,
             range: { start: sweepStart, end: sweepEnd, steps: sweepSteps }
           };
           const response = await ASCENDService.predictBatchSweep(payload);
           if (response.type === 'sweep' && Array.isArray(response.data)) {
             setSweepResults(response.data);
           }
           setSweepIc50Result(null);
           setSweepIc50Error(null);
        } else {
           // Upload mode logic (handled by file input usually, but here triggering 'process' simulation)
           // In a real app, this would submit the uploaded file.
           // We'll assume a dummy file was "ready".
           //const dummyFile = new File([""], "batch.csv");
           //const response = await ASCENDService.processBatchUpload(dummyFile);
// 修改：上面两行注释掉 新增四行 上传真实文件 不再上传空文�?
           if (!selectedFile) {
             throw new Error("请先选择要上传的文件");
           }
           const response = await ASCENDService.processBatchUpload(selectedFile, batchCellLine);
           if (response.type === 'upload' && !Array.isArray(response.data)) {
             setUploadResult(response.data);
           }
        }
      }
      setAppState(AppState.SUCCESS);
    } catch (err: any) {
      console.error(err);
      setErrorMsg(err.message || "An unexpected error occurred.");
      setAppState(AppState.ERROR);
    }
  };

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    // Just visual feedback for the demo
    //console.log("File selected:", e.target.files?.[0]);
//修改:注释上面一行，新增五行
    const file = e.target.files?.[0];//从文件列表里选第一个上�?
    if (file) {
      setSelectedFile(file);//将文件对象存入React状态，动态更�?
      console.log("File selected:", file.name, file.size, "bytes");
    }
  };

  const handleNavigate = (page: 'transcriptomics' | 'metabolomics' | 'enrichment') => {
    const hash = page === 'metabolomics' ? '#metabolomics' : page === 'enrichment' ? '#enrichment' : '#transcriptomics';
    if (window.location.hash !== hash) {
      window.location.hash = hash;
    }
    setActivePage(page);
  };

  const buildTranscriptomicsCsv = (results: GeneExpressionResult[]) => {
    const header = ['cell_id', ...results.map((r) => r.geneId)].join(',');
    const values = results.map((r) => {
      const val = Number.isFinite(r.expressionLevel) ? r.expressionLevel : 0;
      return val.toString();
    });
    const rowId = `sample_${selectedCellLine}_${time}h_${dose}uM`;
    const row = [rowId, ...values].join(',');
    return `${header}\n${row}\n`;
  };

  const buildTranscriptomicsFile = () => {
    if (!singleResults?.length) return null;
    const csv = buildTranscriptomicsCsv(singleResults);
    const safeCell = selectedCellLine.replace(/[^a-zA-Z0-9_-]+/g, '_');
    const filename = `transcriptomics_${safeCell}_${time}h_${dose}uM.csv`;
    const file = new File([csv], filename, { type: 'text/csv' });
    return { file, filename };
  };

  const downloadFile = (file: File) => {
    const url = URL.createObjectURL(file);
    const a = document.createElement('a');
    a.href = url;
    a.download = file.name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const handleExportSweepTranscriptomics = async (action: 'download' | 'send') => {
    if (!sweepResults?.length) {
      setErrorMsg('Sweep results are not ready yet.');
      return;
    }
    const payload: BatchSweepRequest = {
      smiles: batchSmiles,
      cellLineId: batchCellLine,
      sweepVariable: sweepVar,
      fixedParamValue: batchFixedParam,
      range: { start: sweepStart, end: sweepEnd, steps: sweepSteps }
    };
    try {
      const blob = await ASCENDService.fetchSweepMetabolomicsCsv(payload);
      const filename = `metabolomics_sweep_${sweepVar}_${sweepStart}-${sweepEnd}_${sweepSteps}.csv`;
      const file = new File([blob], filename, { type: 'text/csv' });
      setMetabolomicsSeedFile(file);
      setMetabolomicsSeedLabel(filename);
      if (action === 'download') {
        downloadFile(file);
      } else {
        handleNavigate('metabolomics');
      }
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Failed to prepare metabolomics sweep input.');
    }
  };
  const handleExportBatchTranscriptomics = async (action: 'download' | 'send') => {
    if (!uploadResult?.fileId) {
      setErrorMsg('Batch transcriptomics results are not ready yet.');
      return;
    }
    try {
      const blob = await ASCENDService.fetchBatchMetabolomicsCsv(uploadResult.fileId);
      const filename = `metabolomics_${uploadResult.fileId}.csv`;
      const file = new File([blob], filename, { type: 'text/csv' });
      setMetabolomicsSeedFile(file);
      setMetabolomicsSeedLabel(filename);
      if (action === 'download') {
        downloadFile(file);
      } else {
        handleNavigate('metabolomics');
      }
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Failed to prepare metabolomics input.');
    }
  };
  const handleExportTranscriptomics = (action: 'download' | 'send') => {
    const payload = buildTranscriptomicsFile();
    if (!payload) {
      setErrorMsg('Transcriptomics results are not ready yet.');
      return;
    }
    setMetabolomicsSeedFile(payload.file);
    setMetabolomicsSeedLabel(payload.filename);
    if (action === 'download') {
      downloadFile(payload.file);
    } else {
      handleNavigate('metabolomics');
    }
  };

  const handleSendSingleToEnrichment = () => {
    if (!singleResults?.length) {
      setErrorMsg('Transcriptomics results are not ready yet.');
      return;
    }
    setEnrichmentSeedResults(singleResults);
    setEnrichmentSeedLabel(selectedCellLine + ' | ' + time + 'h | ' + dose + 'uM');
    handleNavigate('enrichment');
  };

  const handleSendBatchToEnrichment = async () => {
    if (!uploadResult?.fileId) {
      setErrorMsg('Batch transcriptomics results are not ready yet.');
      return;
    }
    setEnrichmentLoading(true);
    setErrorMsg(null);
    try {
      const data = await ASCENDService.fetchBatchEnrichmentSummary(uploadResult.fileId);
      setEnrichmentSeedResults(data);
      setEnrichmentSeedLabel('Batch Upload | ' + uploadResult.fileId + ' | ' + uploadResult.processedRows + ' samples');
      handleNavigate('enrichment');
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Failed to prepare enrichment input.');
    } finally {
      setEnrichmentLoading(false);
    }
  };

  const handleSendSweepToEnrichment = async () => {
    if (!sweepResults?.length) {
      setErrorMsg('Sweep results are not ready yet.');
      return;
    }
    setEnrichmentLoading(true);
    setErrorMsg(null);
    try {
      const payload: BatchSweepRequest = {
        smiles: batchSmiles,
        cellLineId: batchCellLine,
        sweepVariable: sweepVar,
        fixedParamValue: batchFixedParam,
        range: { start: sweepStart, end: sweepEnd, steps: sweepSteps }
      };
      const data = await ASCENDService.fetchSweepEnrichmentSummary(payload);
      setEnrichmentSeedResults(data);
      setEnrichmentSeedLabel('Sweep ' + sweepVar + ' | ' + sweepStart + '-' + sweepEnd + ' (' + sweepSteps + ' steps)');
      handleNavigate('enrichment');
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Failed to prepare enrichment input.');
    } finally {
      setEnrichmentLoading(false);
    }
  };

  // --- UI Components ---
  const Header = () => (
    <header className="sticky top-0 z-50 w-full border-b border-slate-200 dark:border-slate-800 bg-white/90 dark:bg-[#020617]/90 backdrop-blur-md transition-colors duration-300">
      <div className="container mx-auto px-4 h-16 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <div className="p-2 bg-cyan-50 dark:bg-cyan-950/40 rounded-lg border border-cyan-100 dark:border-cyan-500/20">
            <Dna className="w-6 h-6 text-cyan-600 dark:text-cyan-400" />
          </div>
          <span className="text-xl font-bold tracking-tight text-slate-900 dark:text-white">
            ASCEND <span className="text-slate-400 dark:text-slate-600 text-sm font-mono font-normal ml-2 hidden sm:inline-block">| Virtual Cell AI</span>
          </span>
        </div>
        <div className="flex-1 flex justify-center">
          <div className="flex items-center gap-1 rounded-xl border border-slate-200 dark:border-slate-800 bg-white/80 dark:bg-slate-900/80 p-1">
            <button
              type="button"
              onClick={() => handleNavigate('transcriptomics')}
              className={activePage === 'transcriptomics'
                ? 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all bg-cyan-50 dark:bg-cyan-950/40 text-cyan-700 dark:text-cyan-300 shadow-sm'
                : 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all text-slate-500 hover:text-slate-700 dark:text-slate-400 dark:hover:text-slate-200'}
            >
              Transcriptomics
            </button>
            <button
              type="button"
              onClick={() => handleNavigate('enrichment')}
              className={activePage === 'enrichment'
                ? 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all bg-cyan-50 dark:bg-cyan-950/40 text-cyan-700 dark:text-cyan-300 shadow-sm'
                : 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all text-slate-500 hover:text-slate-700 dark:text-slate-400 dark:hover:text-slate-200'}
            >
              Enrichment
            </button>
            <button
              type="button"
              onClick={() => handleNavigate('metabolomics')}
              className={activePage === 'metabolomics'
                ? 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all bg-cyan-50 dark:bg-cyan-950/40 text-cyan-700 dark:text-cyan-300 shadow-sm'
                : 'px-4 py-1.5 text-xs font-semibold rounded-lg transition-all text-slate-500 hover:text-slate-700 dark:text-slate-400 dark:hover:text-slate-200'}
            >
              Metabolomics
            </button>
          </div>
        </div>
        <div className="flex items-center gap-4">
          <div className="hidden md:flex items-center gap-4 text-xs font-mono text-slate-500 dark:text-slate-400">
             <div className="flex items-center gap-1">
              <span className="w-2 h-2 rounded-full bg-green-500 animate-pulse"></span>
              SYSTEM ONLINE
            </div>
            <span className="px-2 py-1 bg-slate-100 dark:bg-slate-900 rounded border border-slate-200 dark:border-slate-800 text-slate-600 dark:text-slate-400 font-semibold">{MODEL_VERSION}</span>
          </div>
          <button 
            onClick={() => setIsDark(!isDark)}
            className="p-2 rounded-full bg-slate-100 dark:bg-slate-800 text-slate-600 dark:text-slate-400 hover:bg-slate-200 dark:hover:bg-slate-700 transition-colors border border-slate-200 dark:border-slate-700"
          >
            {isDark ? <Sun className="w-5 h-5" /> : <Moon className="w-5 h-5" />}
          </button>
        </div>
      </div>
    </header>
  );

  return (
    <div className="min-h-screen bg-slate-50 dark:bg-[#020617] transition-colors duration-300">
      <Header />

      <main className="container mx-auto px-4 py-8">
        {activePage === 'metabolomics' ? (
          <MetabolomicsPage isDark={isDark} transcriptomicsSeedFile={metabolomicsSeedFile} transcriptomicsSeedLabel={metabolomicsSeedLabel} />
        ) : activePage === 'enrichment' ? (
          <EnrichmentPage isDark={isDark} differentialResults={enrichmentSeedResults || singleResults} contextLabel={enrichmentContextLabel} />
        ) : (
          <>
        
        {/* Intro Section */}
        <div className="mb-10 text-center max-w-3xl mx-auto">
          <h1 className="text-4xl md:text-5xl font-bold text-slate-900 dark:text-white mb-6 tracking-tight transition-colors">
            Predict <span className="text-transparent bg-clip-text bg-gradient-to-r from-cyan-600 to-blue-600 dark:from-cyan-400 dark:to-blue-500">Transcriptomic Perturbations</span>
          </h1>
        </div>

        {/* Mode Switcher */}
        <div className="flex justify-center mb-8">
          <div className="bg-white dark:bg-slate-900 p-1.5 rounded-xl border border-slate-200 dark:border-slate-800 inline-flex items-center gap-1 shadow-sm">
            <button 
              onClick={() => setMode('single')}
              className={`px-6 py-2 rounded-lg text-sm font-medium transition-all flex items-center gap-2 ${mode === 'single' ? 'bg-cyan-50 dark:bg-cyan-950/40 text-cyan-700 dark:text-cyan-400 shadow-sm' : 'text-slate-500 hover:text-slate-700 dark:text-slate-400 dark:hover:text-slate-200'}`}
            >
              <Microscope className="w-4 h-4" />
              Single Analysis
            </button>
            <button 
              onClick={() => setMode('batch')}
              className={`px-6 py-2 rounded-lg text-sm font-medium transition-all flex items-center gap-2 ${mode === 'batch' ? 'bg-cyan-50 dark:bg-cyan-950/40 text-cyan-700 dark:text-cyan-400 shadow-sm' : 'text-slate-500 hover:text-slate-700 dark:text-slate-400 dark:hover:text-slate-200'}`}
            >
              <Layers className="w-4 h-4" />
              Batch Processing
            </button>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-12 gap-8">
          
          {/* LEFT PANEL: CONFIGURATION */}
          <div className="lg:col-span-4 space-y-6">
            <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-xl shadow-slate-200/40 dark:shadow-none relative overflow-hidden transition-all duration-300">
              <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-cyan-500 to-blue-600"></div>
              
              <h2 className="text-lg font-bold text-slate-900 dark:text-white mb-6 flex items-center gap-2">
                <Cpu className="w-5 h-5 text-cyan-600 dark:text-cyan-400" />
                {mode === 'single' ? 'Simulation Parameters' : 'Batch Configuration'}
              </h2>

              <form onSubmit={handlePredict}>
                {mode === 'single' ? (
                  // --- SINGLE MODE FORM ---
                  <>
                    <div className="mb-6">
                      <Input 
                        label="Compound (SMILES)"
                        placeholder="Enter SMILES string..."
                        value={smiles}
                        onChange={(e) => setSmiles(e.target.value)}
                        required
                      />
                      <div className="text-[10px] text-slate-500 mt-1 font-mono break-all bg-slate-50 dark:bg-slate-950/50 p-2 rounded border border-slate-200 dark:border-slate-800">
                        Current: {smiles.substring(0, 30)}...
                      </div>
                    </div>

                    <div className="space-y-6 mb-8 border-t border-slate-100 dark:border-slate-800 pt-6">
                      <RangeSlider 
                        label="Exposure Time"
                        value={time}
                        min={0}
                        max={MAX_TIME_HOURS}
                        unit="h"
                        onChange={(e) => setTime(Number(e.target.value))}
                      />
                      <RangeSlider 
                        label="Dosage"
                        value={dose}
                        min={0}
                        max={MAX_DOSE_UM}
                        unit="μM"
                        onChange={(e) => setDose(Number(e.target.value))}
                      />
                    </div>
                    <div className="mb-8 border-t border-slate-100 dark:border-slate-800 pt-6">
                      <label className="block mb-3 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">
                        Cell Line Origin
                      </label>
                      <div className="flex gap-2 mb-4">
                        <button
                          type="button"
                          onClick={() => setCellLineMode('preset')}
                          className={`flex-1 py-2 text-sm rounded-lg transition-all border font-medium ${cellLineMode === 'preset' ? 'bg-cyan-50 dark:bg-cyan-900/30 text-cyan-700 dark:text-cyan-300 border-cyan-200 dark:border-cyan-700' : 'bg-slate-50 dark:bg-slate-800/50 text-slate-500 dark:text-slate-400 border-slate-200 dark:border-slate-700'}`}
                        >
                          Preset DB
                        </button>
                        <button
                          type="button"
                          onClick={() => setCellLineMode('custom')}
                          className={`flex-1 py-2 text-sm rounded-lg transition-all border font-medium ${cellLineMode === 'custom' ? 'bg-cyan-50 dark:bg-cyan-900/30 text-cyan-700 dark:text-cyan-300 border-cyan-200 dark:border-cyan-700' : 'bg-slate-50 dark:bg-slate-800/50 text-slate-500 dark:text-slate-400 border-slate-200 dark:border-slate-700'}`}
                        >
                          Upload
                        </button>
                      </div>

                      {cellLineMode === 'preset' ? (
                         <div className="relative">
                           <select 
                             value={selectedCellLine}
                             onChange={(e) => setSelectedCellLine(e.target.value)}
                             className="w-full appearance-none bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg px-4 py-3 text-sm text-slate-900 dark:text-white outline-none"
                           >
                             {PRESET_CELL_LINES.map(cell => (
                               <option key={cell.id} value={cell.id}>{cell.name} ({cell.tissue})</option>
                             ))}
                           </select>
                           <ChevronDown className="absolute right-4 top-3.5 w-4 h-4 text-slate-500 pointer-events-none" />
                         </div>
                      ) : (
                        <div className="border-2 border-dashed border-slate-200 dark:border-slate-700 rounded-lg p-4 text-center cursor-pointer hover:border-cyan-500 transition-colors bg-slate-50 dark:bg-slate-950/50">
                          <p className="text-xs text-slate-500">Feature not available in demo</p>
                        </div>
                      )}
                    </div>
                  </>
                ) : (
                  // --- BATCH MODE FORM ---
                  <>
                     <div className="grid grid-cols-2 gap-2 mb-6">
                        <button
                           type="button"
                           onClick={() => setBatchMethod('sweep')}
                           className={`p-3 rounded-lg text-xs font-bold border text-center transition-all ${batchMethod === 'sweep' ? 'border-cyan-500 bg-cyan-50 dark:bg-cyan-900/20 text-cyan-700 dark:text-cyan-400' : 'border-slate-200 dark:border-slate-700 bg-slate-50 dark:bg-slate-800/50 text-slate-500'}`}
                        >
                           Continuous Sweep
                        </button>
                        <button
                           type="button"
                           onClick={() => setBatchMethod('upload')}
                           className={`p-3 rounded-lg text-xs font-bold border text-center transition-all ${batchMethod === 'upload' ? 'border-cyan-500 bg-cyan-50 dark:bg-cyan-900/20 text-cyan-700 dark:text-cyan-400' : 'border-slate-200 dark:border-slate-700 bg-slate-50 dark:bg-slate-800/50 text-slate-500'}`}
                        >
                           Batch Upload
                        </button>
                     </div>

                     {batchMethod === 'sweep' ? (
                       <div className="space-y-6">
                          <div>
                            <label className="block mb-2 text-xs font-mono font-bold text-slate-500 dark:text-cyan-500 uppercase">Variable Parameter</label>
                            <div className="flex gap-2">
                               <button type="button" onClick={() => setSweepVar('dose')} className={`flex-1 py-2 text-xs rounded border ${sweepVar === 'dose' ? 'bg-blue-500 text-white border-blue-600' : 'bg-slate-100 dark:bg-slate-800 border-slate-200 dark:border-slate-700 text-slate-500'}`}>Dose (μM)</button>
                               <button type="button" onClick={() => setSweepVar('time')} className={`flex-1 py-2 text-xs rounded border ${sweepVar === 'time' ? 'bg-blue-500 text-white border-blue-600' : 'bg-slate-100 dark:bg-slate-800 border-slate-200 dark:border-slate-700 text-slate-500'}`}>Time (h)</button>
                            </div>
                          </div>

                          <div className="grid grid-cols-3 gap-2">
                             <Input label="Start" type="number" value={sweepStart} onChange={e => setSweepStart(Number(e.target.value))} />
                             <Input label="End" type="number" value={sweepEnd} onChange={e => setSweepEnd(Number(e.target.value))} />
                             <Input label="Steps" type="number" value={sweepSteps} max={20} onChange={e => setSweepSteps(Number(e.target.value))} />
                          </div>

                          <div className="pt-4 border-t border-slate-200 dark:border-slate-800">
                             <p className="text-xs font-bold text-slate-400 mb-3 uppercase">Fixed Parameters</p>
                             <Input 
                                label={sweepVar === 'dose' ? 'Fixed Time (h)' : 'Fixed Dose (μM)'} 
                                type="number" 
                                value={batchFixedParam} 
                                onChange={e => setBatchFixedParam(Number(e.target.value))} 
                             />
                             <Input 
                                label="Compound" 
                                value={batchSmiles} 
                                onChange={e => setBatchSmiles(e.target.value)} 
                             />
                             <div className="relative mt-4">
                               <select 
                                 value={batchCellLine}
                                 onChange={(e) => setBatchCellLine(e.target.value)}
                                 className="w-full appearance-none bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg px-4 py-3 text-sm text-slate-900 dark:text-white outline-none"
                               >
                                 {PRESET_CELL_LINES.map(cell => (
                                   <option key={cell.id} value={cell.id}>{cell.name}</option>
                                 ))}
                               </select>
                               <ChevronDown className="absolute right-4 top-3.5 w-4 h-4 text-slate-500 pointer-events-none" />
                             </div>
                          </div>
                       </div>
                     ) : (
                       <div className="space-y-4">
                          <div className="p-4 bg-slate-50 dark:bg-slate-950/30 rounded-lg border border-slate-200 dark:border-slate-800 text-sm text-slate-600 dark:text-slate-400">
                             <p className="mb-2">Upload a CSV file containing multiple queries. Requires columns: <code className="bg-slate-200 dark:bg-slate-800 px-1 rounded">smiles</code>, <code className="bg-slate-200 dark:bg-slate-800 px-1 rounded">dose</code>, <code className="bg-slate-200 dark:bg-slate-800 px-1 rounded">time</code>.</p>
                             <button 
                                type="button" 
                                onClick={() => {
                                  const templateContent = "smiles,time,dose\nCCO,24,10\nCC1=CC=CC=C1,12,5\nCCN(CC)CC,48,20";
                                  const blob = new Blob([templateContent], { type: 'text/csv' });
                                  const url = window.URL.createObjectURL(blob);
                                  const a = document.createElement('a');
                                  a.href = url;
                                  a.download = 'batch_template.csv';
                                  document.body.appendChild(a);
                                  a.click();
                                  document.body.removeChild(a);
                                  window.URL.revokeObjectURL(url);
                                }}
                                className="text-cyan-600 hover:underline flex items-center gap-1 text-xs font-bold mt-2"
                             >
                                <Download className="w-3 h-3" /> Download Template
                             </button>
                          </div>

                          <div className="border-2 border-dashed border-slate-200 dark:border-slate-700 rounded-xl p-8 text-center hover:border-cyan-400 dark:hover:border-cyan-500 transition-colors bg-white dark:bg-slate-950/50 cursor-pointer group relative">
                            {selectedFile ? (
                              <>
                                <FileSpreadsheet className="w-10 h-10 text-green-500 transition-colors mx-auto mb-3" />
                                <p className="text-sm text-slate-700 dark:text-slate-300 mb-2 font-medium">{selectedFile.name}</p>
                                <p className="text-xs text-slate-500">{(selectedFile.size / 1024).toFixed(1)} KB</p>
                                <button 
                                  type="button" 
                                  onClick={(e) => {
                                    e.stopPropagation();
                                    setSelectedFile(null);
                                  }}
                                  className="mt-2 text-xs text-red-500 hover:text-red-700 font-medium"
                                >
                                  Remove File
                                </button>
                              </>
                            ) : (
                              <>
                                <FileSpreadsheet className="w-10 h-10 text-slate-300 group-hover:text-cyan-500 transition-colors mx-auto mb-3" />
                                <p className="text-sm text-slate-500 dark:text-slate-400 mb-2 font-medium">Drag & Drop or Click to Upload</p>
                                <p className="text-xs text-slate-400">.CSV, .XLSX (Max 50MB)</p>
                              </>
                            )}
                            <input type="file" className="hidden" id="batch-upload" onChange={handleFileUpload} />
                            <label htmlFor="batch-upload" className="absolute inset-0 cursor-pointer"></label>
                          </div>
                          <div className="pt-2">
                            <label className="block mb-2 text-xs font-mono font-bold text-slate-500 dark:text-cyan-500 uppercase">
                              Cell Line
                            </label>
                            <div className="relative">
                              <select
                                value={batchCellLine}
                                onChange={(e) => setBatchCellLine(e.target.value)}
                                className="w-full appearance-none bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg px-4 py-3 text-sm text-slate-900 dark:text-white outline-none"
                              >
                                {PRESET_CELL_LINES.map(cell => (
                                  <option key={cell.id} value={cell.id}>{cell.name}</option>
                                ))}
                              </select>
                              <ChevronDown className="absolute right-4 top-3.5 w-4 h-4 text-slate-500 pointer-events-none" />
                            </div>
                          </div>
                       </div>
                     )}
                  </>
                )}

                <Button 
                  type="submit" 
                  className="w-full mt-8 shadow-lg shadow-cyan-500/20 hover:shadow-cyan-500/40" 
                  isLoading={appState === AppState.LOADING}
                >
                  {appState === AppState.LOADING ? 'Processing...' : (mode === 'single' ? 'Generate Prediction' : (batchMethod === 'sweep' ? 'Run Parameter Sweep' : 'Process Batch File'))}
                </Button>

                {errorMsg && (
                  <div className="mt-4 p-3 bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/30 rounded-lg text-red-600 dark:text-red-400 text-xs">
                    {errorMsg}
                  </div>
                )}
              </form>
            </div>
          </div>

          {/* RIGHT PANEL: VISUALIZATION */}
          <div className="lg:col-span-8 space-y-6 min-w-0">
            
            {appState === AppState.IDLE && (
              <div className="h-full min-h-[500px] flex flex-col items-center justify-center bg-white dark:bg-slate-900 border border-dashed border-slate-300 dark:border-slate-800 rounded-xl transition-colors shadow-sm">
                <div className="p-4 bg-slate-50 dark:bg-slate-800 rounded-full mb-4">
                  <Microscope className="w-12 h-12 text-slate-400 dark:text-slate-500" />
                </div>
                <h3 className="text-slate-900 dark:text-white font-medium text-lg">Ready for Simulation</h3>
                <p className="text-slate-500 dark:text-slate-400 text-sm mt-2 max-w-md text-center">
                  {mode === 'single' ? 'Configure parameters on the left to initialize the ASCEND Virtual Cell.' : 'Setup a parameter sweep or upload a batch file to begin processing multiple simulations.'}
                </p>
              </div>
            )}

            {appState === AppState.LOADING && (
              <div className="h-full min-h-[500px] flex flex-col items-center justify-center bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl relative overflow-hidden transition-colors shadow-lg">
                <div className="absolute top-0 w-full h-1 bg-cyan-500/50 shadow-[0_0_15px_rgba(6,182,212,0.5)] animate-scan"></div>
                <div className="relative">
                   <div className="absolute inset-0 bg-cyan-400/20 blur-xl rounded-full"></div>
                   <Activity className="relative w-16 h-16 text-cyan-500 animate-pulse mb-6" />
                </div>
                <h3 className="text-slate-900 dark:text-white font-mono text-lg font-bold">PROCESSING DATA</h3>
                <p className="text-slate-500 text-xs mt-2 font-mono">{mode === 'batch' ? 'Running batch inference matrix...' : 'Calculating transcriptomic perturbations...'}</p>
              </div>
            )}

            {appState === AppState.SUCCESS && (
               <div className="animate-fade-in space-y-6">
                 
                 {mode === 'single' && singleResults && (
                   <>
                    {/* Single Mode Stats */}
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                      <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-5 rounded-xl flex items-center gap-4 shadow-md shadow-slate-200/50 dark:shadow-none transition-all">
                        <div className="p-3 bg-purple-50 dark:bg-purple-900/20 rounded-lg text-purple-600 dark:text-purple-400"><FlaskConical className="w-6 h-6"/></div>
                        <div>
                          <div className="text-xs text-slate-500 font-bold uppercase font-mono tracking-wider">Input Dose</div>
                          <div className="text-xl font-bold text-slate-900 dark:text-white">{dose} μM</div>
                        </div>
                      </div>
                      <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-5 rounded-xl flex items-center gap-4 shadow-md shadow-slate-200/50 dark:shadow-none transition-all">
                        <div className="p-3 bg-blue-50 dark:bg-blue-900/20 rounded-lg text-blue-600 dark:text-blue-400"><Clock className="w-6 h-6"/></div>
                        <div>
                          <div className="text-xs text-slate-500 font-bold uppercase font-mono tracking-wider">Duration</div>
                          <div className="text-xl font-bold text-slate-900 dark:text-white">{time} Hours</div>
                        </div>
                      </div>
                      <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-5 rounded-xl flex items-center gap-4 shadow-md shadow-slate-200/50 dark:shadow-none transition-all">
                        <div className="p-3 bg-green-50 dark:bg-green-900/20 rounded-lg text-green-600 dark:text-green-400"><Zap className="w-6 h-6"/></div>
                        <div>
                          <div className="text-xs text-slate-500 font-bold uppercase font-mono tracking-wider">Cell Line</div>
                          <div className="text-xl font-bold text-slate-900 dark:text-white">{selectedCellLine}</div>
                        </div>
                      </div>
                    </div>
                    {/* Single Mode Chart */}
                    <div className="relative">
                       <ResultsChart data={singleResults} isDark={isDark} />
                    </div>
                    <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6 rounded-xl shadow-md space-y-4">
                      <div className="flex flex-col lg:flex-row lg:items-center lg:justify-between gap-4">
                        <div>
                          <h3 className="text-lg font-bold text-slate-900 dark:text-white">Transcriptomics Output</h3>
                          <p className="text-sm text-slate-500 dark:text-slate-400">
                            Download the full gene expression prediction or send it directly to the Metabolomics module.
                          </p>
                        </div>
                        <div className="flex flex-wrap gap-3">
                          <Button
                            variant="secondary"
                            className="flex items-center gap-2"
                            onClick={() => handleExportTranscriptomics('download')}
                          >
                            <Download className="w-4 h-4" /> Download CSV
                          </Button>
                          <Button
                            variant="outline"
                            className="flex items-center gap-2"
                            onClick={handleSendSingleToEnrichment}
                          >
                            <Dna className="w-4 h-4" /> Send to Enrichment
                          </Button>
                          <Button
                            variant="outline"
                            className="flex items-center gap-2"
                            onClick={() => handleExportTranscriptomics('send')}
                          >
                            <FlaskConical className="w-4 h-4" /> Send to Metabolomics
                          </Button>
                        </div>
                      </div>
                    </div>
                    <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6 rounded-xl shadow-md space-y-4">
                      <div className="flex items-center justify-between gap-4">
                        <div>
                          <h3 className="text-lg font-bold text-slate-900 dark:text-white">
                            Module 2: IC50 Prediction
                          </h3>
                          <p className="text-sm text-slate-500 dark:text-slate-400">
                            Use the transcriptomic signature from Module 1 as input to estimate IC50.
                          </p>
                        </div>
                        <Button onClick={handleIC50Prediction} isLoading={ic50Loading}>
                          {ic50Result ? 'Re-run IC50' : 'Run IC50'}
                        </Button>
                      </div>

                      {ic50Result ? (
                        <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
                          <div className="p-4 rounded-lg bg-gradient-to-r from-cyan-500 to-blue-500 text-white">
                            <div className="flex items-center justify-between gap-3">
                              <p className="text-xs uppercase tracking-widest opacity-80">
                                {ic50Display === 'ic50' ? 'IC50' : 'ln(IC50)'}
                              </p>
                              <div className="flex items-center gap-1 rounded-full bg-white/15 p-1">
                                <button
                                  type="button"
                                  onClick={() => setIc50Display('ic50')}
                                  className={`px-2 py-0.5 text-[10px] rounded-full ${ic50Display === 'ic50' ? 'bg-white text-slate-900' : 'text-white/80 hover:text-white'}`}
                                >
                                  IC50
                                </button>
                                <button
                                  type="button"
                                  onClick={() => setIc50Display('lnIc50')}
                                  className={`px-2 py-0.5 text-[10px] rounded-full ${ic50Display === 'lnIc50' ? 'bg-white text-slate-900' : 'text-white/80 hover:text-white'}`}
                                >
                                  ln(IC50)
                                </button>
                              </div>
                            </div>
                            <p className="text-4xl font-bold">
                              {ic50Display === 'ic50'
                                ? `${ic50Result.ic50.toFixed(4)} ${ic50Result.unit}`
                                : ic50Result.lnIc50.toFixed(4)}
                            </p>
                          </div>
                          {ic50Result.confidence !== undefined && (
                            <div className="p-4 rounded-lg bg-slate-100 dark:bg-slate-800">
                              <p className="text-xs uppercase tracking-widest text-slate-500 dark:text-slate-400">Confidence</p>
                              <p className="text-2xl font-semibold text-slate-900 dark:text-white">
                                {(ic50Result.confidence * 100).toFixed(1)}%
                              </p>
                            </div>
                          )}
                          <div className="p-4 rounded-lg border border-slate-200 dark:border-slate-700">
                            <p className="text-xs uppercase tracking-widest text-slate-500 dark:text-slate-400">Tool</p>
                            <p className="text-sm font-mono text-slate-700 dark:text-slate-300">{ic50Result.toolVersion ?? 'custom'}</p>
                          </div>
                        </div>
                      ) : (
                        <p className="text-sm text-slate-500 dark:text-slate-400">
                          Module 1 completed. Click "Run IC50" to start prediction.
                        </p>
                      )}

                      {ic50Error && (
                        <div className="p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
                          {ic50Error}
                        </div>
                      )}
                    </div>
                   </>
                 )}

                 {mode === 'batch' && batchMethod === 'sweep' && sweepResults && (
                   <>
                      <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6 rounded-xl shadow-md">
                         <h3 className="text-lg font-bold text-slate-900 dark:text-white mb-2">Parameter Sweep Complete</h3>
                         <p className="text-slate-500 text-sm">Visualizing gene expression response across {sweepVar} range {sweepStart}-{sweepEnd}.</p>
                      </div>
                      <div className="relative">
                         <SweepLineChart data={sweepResults} isDark={isDark} xAxisLabel={sweepVar === 'dose' ? 'Concentration (μM)' : 'Time (Hours)'} />
                      </div>
                      <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6 rounded-xl shadow-md space-y-4">
                        <div className="flex flex-col lg:flex-row lg:items-center lg:justify-between gap-4">
                          <div>
                            <h3 className="text-lg font-bold text-slate-900 dark:text-white">Metabolomics Input (Sweep)</h3>
                            <p className="text-sm text-slate-500 dark:text-slate-400">
                              Export the full gene matrix for all sweep points to use in metabolomics prediction.
                            </p>
                          </div>
                          <div className="flex flex-wrap gap-3">

                            <Button
                              variant="outline"
                              className="flex items-center gap-2"
                              onClick={() => handleExportSweepTranscriptomics('send')}
                            >
                              <FlaskConical className="w-4 h-4" /> Send to Metabolomics
                            </Button>
                            <Button
                              variant="outline"
                              className="flex items-center gap-2"
                              onClick={handleSendSweepToEnrichment}
                              isLoading={enrichmentLoading}
                            >
                              <Dna className="w-4 h-4" /> Send to Enrichment
                            </Button>
                          </div>
                        </div>
                      </div>
                      <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6 rounded-xl shadow-md space-y-4">
                        <div className="flex items-center justify-between gap-4">
                          <div>
                            <h3 className="text-lg font-bold text-slate-900 dark:text-white">
                              Module 2: Sweep IC50 Prediction
                            </h3>
                            <p className="text-sm text-slate-500 dark:text-slate-400">
                              IC50 is calculated once per unique SMILES + Cell Line pair across sweep points.
                            </p>
                          </div>
                          <Button onClick={handleSweepIC50Prediction} isLoading={sweepIc50Loading}>
                            {sweepIc50Result ? 'Re-run IC50 Sweep' : 'Run IC50 Sweep'}
                          </Button>
                        </div>

                        {sweepIc50Result ? (
                          <>
                            <div className="overflow-x-auto border border-slate-200 dark:border-slate-800 rounded-lg">
                              <table className="min-w-full text-left text-sm">
                                <thead className="bg-slate-50 dark:bg-slate-800 text-slate-500 dark:text-slate-400 uppercase text-xs">
                                  <tr>
                                    <th className="px-4 py-2">#</th>
                                    <th className="px-4 py-2">SMILES</th>
                                    <th className="px-4 py-2">Time (h)</th>
                                    <th className="px-4 py-2">Dose (uM)</th>
                                    <th className="px-4 py-2">IC50</th>
                                    <th className="px-4 py-2">ln(IC50)</th>
                                    <th className="px-4 py-2">Confidence</th>
                                  </tr>
                                </thead>
                                <tbody className="text-slate-700 dark:text-slate-300">
                                  {sweepIc50Result.preview.map((row, idx) => (
                                    <tr key={`${row.smiles}-${idx}`} className="border-t border-slate-200 dark:border-slate-800">
                                      <td className="px-4 py-2">{idx + 1}</td>
                                      <td className="px-4 py-2 max-w-[240px] truncate" title={row.smiles}>{row.smiles}</td>
                                      <td className="px-4 py-2">{row.time}</td>
                                      <td className="px-4 py-2">{row.dose}</td>
                                      <td className="px-4 py-2">{row.ic50.toFixed(4)} {row.unit}</td>
                                      <td className="px-4 py-2">{row.lnIc50.toFixed(4)}</td>
                                      <td className="px-4 py-2">{row.confidence !== undefined ? `${(row.confidence * 100).toFixed(1)}%` : 'N/A'}</td>
                                    </tr>
                                  ))}
                                </tbody>
                              </table>
                            </div>

                            <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-3">
                              <p className="text-xs text-slate-500 dark:text-slate-400">
                                Showing {sweepIc50Result.preview.length} of {sweepIc50Result.processedRows} results.
                              </p>
                              <Button
                                className="flex items-center gap-2"
                                onClick={() => {
                                  if (sweepIc50Result?.downloadUrl) {
                                    const a = document.createElement('a');
                                    a.href = sweepIc50Result.downloadUrl;
                                    a.download = 'ic50_sweep_results.zip';
                                    document.body.appendChild(a);
                                    a.click();
                                    document.body.removeChild(a);
                                  }
                                }}
                              >
                                <Download className="w-4 h-4" /> Download IC50 Results (.ZIP)
                              </Button>
                            </div>
                          </>
                        ) : (
                          <p className="text-sm text-slate-500 dark:text-slate-400">
                            Sweep completed. Click "Run IC50 Sweep" to generate IC50 predictions.
                          </p>
                        )}

                        {sweepIc50Error && (
                          <div className="p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
                            {sweepIc50Error}
                          </div>
                        )}
                      </div>
                   </>
                 )}

                 {mode === 'batch' && batchMethod === 'upload' && uploadResult && (
                   <>
                     <div className="flex flex-col items-center justify-center min-h-[400px] bg-white dark:bg-slate-900 rounded-xl border border-slate-200 dark:border-slate-800 shadow-md p-10 text-center">
                        <div className="w-20 h-20 bg-green-100 dark:bg-green-900/30 rounded-full flex items-center justify-center mb-6">
                          <CheckCircle2 className="w-10 h-10 text-green-600 dark:text-green-400" />
                        </div>
                        <h2 className="text-2xl font-bold text-slate-900 dark:text-white mb-2">Batch Processing Successful</h2>
                        <p className="text-slate-500 dark:text-slate-400 mb-8 max-w-md">
                          Successfully processed <strong>{uploadResult.processedRows}</strong> samples. The transcriptomic signatures have been generated and compiled.
                        </p>
                        <Button 
                          className="flex items-center gap-2"
                          onClick={() => {
                            if (uploadResult?.downloadUrl) {
                              const a = document.createElement('a');
                              a.href = uploadResult.downloadUrl;
                              a.download = 'batch_results.zip';
                              document.body.appendChild(a);
                              a.click();
                              document.body.removeChild(a);
                            }
                          }}
                        >
                           <Download className="w-4 h-4" /> Download Results (.ZIP)
                        </Button>
                        <div className="mt-6 flex flex-wrap gap-3 justify-center">

                          <Button
                            variant="outline"
                            className="flex items-center gap-2"
                            onClick={() => handleExportBatchTranscriptomics('send')}
                          >
                            <FlaskConical className="w-4 h-4" /> Send to Metabolomics
                          </Button>
                          <Button
                            variant="outline"
                            className="flex items-center gap-2"
                            onClick={handleSendBatchToEnrichment}
                            isLoading={enrichmentLoading}
                          >
                            <Dna className="w-4 h-4" /> Send to Enrichment
                          </Button>
                        </div>
                        {enrichmentLoading && (
                          <p className="text-xs text-slate-500 dark:text-slate-400 mt-2">Preparing enrichment input...</p>
                        )}
                     </div>

                     <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6 rounded-xl shadow-md space-y-4">
                       <div className="flex items-center justify-between gap-4">
                         <div>
                           <h3 className="text-lg font-bold text-slate-900 dark:text-white">
                             Module 2: Batch IC50 Prediction
                           </h3>
                           <p className="text-sm text-slate-500 dark:text-slate-400">
                             Run IC50 inference for all transcriptomic signatures and review the top 5 results.
                           </p>
                         </div>
                         <Button onClick={handleBatchIC50Prediction} isLoading={batchIc50Loading}>
                           {batchIc50Result ? 'Re-run IC50 Batch' : 'Run IC50 Batch'}
                         </Button>
                       </div>

                       {batchIc50Result ? (
                         <>
                           <div className="overflow-x-auto border border-slate-200 dark:border-slate-800 rounded-lg">
                             <table className="min-w-full text-left text-sm">
                               <thead className="bg-slate-50 dark:bg-slate-800 text-slate-500 dark:text-slate-400 uppercase text-xs">
                                 <tr>
                                   <th className="px-4 py-2">#</th>
                                   <th className="px-4 py-2">SMILES</th>
                                   <th className="px-4 py-2">Time (h)</th>
                                   <th className="px-4 py-2">Dose (uM)</th>
                                   <th className="px-4 py-2">IC50</th>
                                   <th className="px-4 py-2">ln(IC50)</th>
                                   <th className="px-4 py-2">Confidence</th>
                                 </tr>
                               </thead>
                               <tbody className="text-slate-700 dark:text-slate-300">
                                 {batchIc50Result.preview.map((row, idx) => (
                                   <tr key={`${row.smiles}-${idx}`} className="border-t border-slate-200 dark:border-slate-800">
                                     <td className="px-4 py-2">{idx + 1}</td>
                                     <td className="px-4 py-2 max-w-[240px] truncate" title={row.smiles}>{row.smiles}</td>
                                     <td className="px-4 py-2">{row.time}</td>
                                     <td className="px-4 py-2">{row.dose}</td>
                                     <td className="px-4 py-2">{row.ic50.toFixed(4)} {row.unit}</td>
                                     <td className="px-4 py-2">{row.lnIc50.toFixed(4)}</td>
                                     <td className="px-4 py-2">{row.confidence !== undefined ? `${(row.confidence * 100).toFixed(1)}%` : 'N/A'}</td>
                                   </tr>
                                 ))}
                               </tbody>
                             </table>
                           </div>

                           <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-3">
                             <p className="text-xs text-slate-500 dark:text-slate-400">
                               Showing {batchIc50Result.preview.length} of {batchIc50Result.processedRows} results.
                             </p>
                             <Button
                               className="flex items-center gap-2"
                               onClick={() => {
                                 if (batchIc50Result?.downloadUrl) {
                                   const a = document.createElement('a');
                                   a.href = batchIc50Result.downloadUrl;
                                   a.download = 'ic50_batch_results.zip';
                                   document.body.appendChild(a);
                                   a.click();
                                   document.body.removeChild(a);
                                 }
                               }}
                             >
                               <Download className="w-4 h-4" /> Download IC50 Results (.ZIP)
                             </Button>
                           </div>
                         </>
                       ) : (
                         <p className="text-sm text-slate-500 dark:text-slate-400">
                           Module 1 completed. Click "Run IC50 Batch" to generate IC50 predictions.
                         </p>
                       )}

                       {batchIc50Error && (
                         <div className="p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
                           {batchIc50Error}
                         </div>
                       )}
                     </div>
                   </>
                 )}

               </div>
            )}

          </div>
        </div>
          </>
        )}
      </main>
    </div>
  );
};

export default App;
