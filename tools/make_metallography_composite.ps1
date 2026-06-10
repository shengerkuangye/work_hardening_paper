Add-Type -AssemblyName System.Drawing

$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
$Src = Join-Path $Root "data\metallography_from_pptx\renamed_metallography_images"
$Out = Join-Path $Root "figures\ta4_metallography_composite.png"

$Samples = @(
    @{ Label = "a"; Stem = "ta4_m" },
    @{ Label = "b"; Stem = "ta4_y_6_48" },
    @{ Label = "c"; Stem = "ta4_y_6_02" },
    @{ Label = "d"; Stem = "ta4_y_5_6" },
    @{ Label = "e"; Stem = "ta4_y_5_25" },
    @{ Label = "f"; Stem = "ta4_y_5" }
)
$Mags = @("100x", "200x", "500x")

$CellW = 430
$CellH = 327
$Gap = 8
$Margin = 8

$Cols = $Samples.Count
$Rows = $Mags.Count
$Width = $Margin * 2 + $Cols * $CellW + ($Cols - 1) * $Gap
$Height = $Margin * 2 + $Rows * $CellH + ($Rows - 1) * $Gap

$Bitmap = New-Object System.Drawing.Bitmap($Width, $Height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$Bitmap.SetResolution(300, 300)
$Graphics = [System.Drawing.Graphics]::FromImage($Bitmap)
$Graphics.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::HighQuality
$Graphics.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
$Graphics.Clear([System.Drawing.Color]::Transparent)

$LabelFont = New-Object System.Drawing.Font("Arial", 28, [System.Drawing.FontStyle]::Bold, [System.Drawing.GraphicsUnit]::Pixel)
$SmallFont = New-Object System.Drawing.Font("Arial", 22, [System.Drawing.FontStyle]::Bold, [System.Drawing.GraphicsUnit]::Pixel)
$Brush = [System.Drawing.Brushes]::Black
$WhiteBrush = [System.Drawing.Brushes]::White
$Pen = New-Object System.Drawing.Pen([System.Drawing.Color]::Black, 2)
$ArrowPen = New-Object System.Drawing.Pen([System.Drawing.Color]::Black, 4)

function Get-CoverSourceRect($Image, [int]$TargetW, [int]$TargetH) {
    $SrcRatio = $Image.Width / $Image.Height
    $DstRatio = $TargetW / $TargetH
    if ($SrcRatio -gt $DstRatio) {
        $NewW = [int]($Image.Height * $DstRatio)
        $Left = [int](($Image.Width - $NewW) / 2)
        return New-Object System.Drawing.Rectangle($Left, 0, $NewW, $Image.Height)
    }
    $NewH = [int]($Image.Width / $DstRatio)
    $Top = [int](($Image.Height - $NewH) / 2)
    return New-Object System.Drawing.Rectangle(0, $Top, $Image.Width, $NewH)
}

for ($r = 0; $r -lt $Rows; $r++) {
    $Y = $Margin + $r * ($CellH + $Gap)

    for ($c = 0; $c -lt $Cols; $c++) {
        $X = $Margin + $c * ($CellW + $Gap)
        $Path = Join-Path $Src "$($Samples[$c].Stem)_z_$($Mags[$r]).jpeg"
        $Image = [System.Drawing.Image]::FromFile($Path)
        $SrcRect = Get-CoverSourceRect $Image $CellW $CellH
        $DstRect = New-Object System.Drawing.Rectangle($X, $Y, $CellW, $CellH)
        $Graphics.DrawImage($Image, $DstRect, $SrcRect, [System.Drawing.GraphicsUnit]::Pixel)
        $Graphics.DrawRectangle($Pen, $X, $Y, $CellW, $CellH)
        $Image.Dispose()

        $SubLabel = "$($Samples[$c].Label)-$($r + 1)"
        $LabelSize = $Graphics.MeasureString($SubLabel, $LabelFont)
        $LabelRect = New-Object System.Drawing.RectangleF(($X + 8), ($Y + 8), ($LabelSize.Width + 16), ($LabelSize.Height + 8))
        $Graphics.FillRectangle($WhiteBrush, $LabelRect)
        $Graphics.DrawRectangle($Pen, [int]$LabelRect.X, [int]$LabelRect.Y, [int]$LabelRect.Width, [int]$LabelRect.Height)
        $Graphics.DrawString($SubLabel, $LabelFont, $Brush, [System.Drawing.PointF]::new(($X + 16), ($Y + 10)))

        # RD marks the deformation direction in the longitudinal section.
        $ArrowX = $X + $CellW - 30
        $ArrowY0 = $Y + $CellH - 98
        $ArrowY1 = $Y + $CellH - 34
        $Graphics.DrawLine($ArrowPen, $ArrowX, $ArrowY1, $ArrowX, $ArrowY0)
        $ArrowHead = @(
            [System.Drawing.Point]::new($ArrowX, $ArrowY0 - 16),
            [System.Drawing.Point]::new($ArrowX - 12, $ArrowY0 + 7),
            [System.Drawing.Point]::new($ArrowX + 12, $ArrowY0 + 7)
        )
        $Graphics.FillPolygon($Brush, $ArrowHead)
        $RdSize = $Graphics.MeasureString("RD", $SmallFont)
        $RdRect = New-Object System.Drawing.RectangleF(($ArrowX - 18), ($ArrowY1 + 4), ($RdSize.Width + 8), ($RdSize.Height + 4))
        $Graphics.FillRectangle($WhiteBrush, $RdRect)
        $Graphics.DrawString("RD", $SmallFont, $Brush, [System.Drawing.PointF]::new(($ArrowX - 14), ($ArrowY1 + 4)))
    }
}

$Bitmap.Save($Out, [System.Drawing.Imaging.ImageFormat]::Png)
$Graphics.Dispose()
$Bitmap.Dispose()

Write-Output $Out
