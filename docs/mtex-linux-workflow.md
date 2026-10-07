# Linux MTEX 受控运行工作流

## 目的

避免 MATLAB/MTEX 已完成数据处理但主进程、图形进程或子进程继续占用 CPU/内存，也避免多个分析批次同时写入同一结果目录。该工作流只清理由本次命令创建的 Linux 会话，不终止预先存在的 MATLAB 或其他用户进程。

## 环境基线

- MATLAB：`/usr/local/MATLAB/R2024a/bin/matlab`
- MTEX：`/home/abcd/MATLAB/toolboxes/mtex-6.1.1`
- 用户启动脚本：`/home/abcd/Documents/MATLAB/startup.m`，会自动加载 MTEX。
- 项目脚本：`tools/mtex/`

不要在批处理表达式中重复调用 `startup_mtex`。不要使用 `matlab ... &`、`nohup matlab ...` 或脱离项目锁的后台会话。

## 标准调用

```bash
tools/run_mtex_guarded.sh --timeout 7200 \
  "run_comprehensive_ebsd_analysis(pwd,fullfile(pwd,'results','mtex_ebsd_comprehensive'))"
```

针对单个生成器：

```bash
tools/run_mtex_guarded.sh --timeout 3600 \
  "generate_c_axis_pole_figures(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans'),fullfile(pwd,'results','mtex_c_axis_pf'))"
```

环境变量可覆盖默认值：

```bash
MATLAB_BIN=/path/to/matlab \
MTEX_TIMEOUT_SECONDS=7200 \
MTEX_COMPLETION_GRACE_SECONDS=20 \
tools/run_mtex_guarded.sh "MATLAB expression"
```

## 运行保证

启动器执行以下操作：

1. 使用 `flock` 获取项目级互斥锁，拒绝重叠运行。
2. 使用 `setsid` 为本次 MATLAB 建立唯一 Linux 会话并记录 PID、PGID 和 SID。
3. 在 `.codex_tmp/mtex_runs/<run-id>/` 保存表达式、stdout、stderr 和运行前后进程清单。
4. MATLAB 表达式正常返回后，从 `MTEX_RUN_DIR` 环境变量重新取得运行目录，再写入 `complete.marker` 并调用 `exit(0)`；即使用户脚本执行了 `clear`，也不会丢失收尾路径。
5. 若完成标记出现后 MATLAB 在宽限期内仍未退出，仅终止同一 SID 下的本次运行进程。
6. 若超过硬超时仍未完成，终止本次 SID，返回退出码 124。
7. 运行结束后再次枚举本次 SID；若仍有进程，返回退出码 70。

因此，清理目标由本次运行的 SID 确定，不根据模糊的进程名全局执行 `pkill` 或 `killall`。

系统中可能预先存在独立 SID 的全局 `MathWorksServiceHost service`。它不等同于未退出的 MTEX 批次，也不属于受控运行的清理范围；除非另有经过进程归属核验的管理任务，不得由本工作流终止。

## 结果与临时文件规定

- 正式派生结果写入 `results/<analysis-name>/`。
- 测试、日志、标记和进程快照写入 `.codex_tmp/`，不得提交 Git。
- 不覆盖 `data/` 下的 CTF、拉伸原始文件或金相源文件。
- 每个正式结果目录应包含可追溯的 CSV/参数文件；PNG/PDF/TIF 不能替代数值表。
- 失败后先读取该次运行目录内的 `stderr.log` 和 `run_metadata.txt`，不要盲目重新运行。

## 返回码

| 返回码 | 含义 |
|---:|---|
| 0 | 表达式完成、完成标记存在、运行 SID 无残留进程 |
| 1 或 MATLAB 返回码 | MATLAB 表达式失败 |
| 70 | 清理后本次 SID 仍有残留进程 |
| 75 | 已有另一个受控 MTEX 批次持有项目锁 |
| 124 | 达到硬超时，已清理本次运行会话 |
| 127 | MATLAB、`flock` 或 `setsid` 不可用 |
